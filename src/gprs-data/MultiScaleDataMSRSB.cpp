#include "MultiScaleDataMSRSB.hpp"
#include "MetisInterface.hpp"
#include "angem/CollisionGJK.hpp"  // collision_gjk algorithm
#include "angem/Collisions.hpp"  // point_inside_surface
#include "mesh/SurfaceMesh.hpp"  // to store support bounding surface
#include "VTKWriter.hpp" // debug bounding region

#include <unordered_set>

namespace multiscale {

using Point = angem::Point<3,double>;
using std::unordered_set;

MultiScaleDataMSRSB::MultiScaleDataMSRSB(mesh::Mesh  & grid,
                                         const size_t  n_blocks)
    :
    grid(grid),
    active_layer_index(0)
{
  auto & layer = layers.emplace_back();
  layer.index = 0;
  layer.n_blocks = n_blocks;
  layer.n_cells = grid.n_cells();
}


void MultiScaleDataMSRSB::build_data()
{
  std::cout << "building METIS partitioning...";
  build_partitioning();
  std::cout << "OK" << std::endl;

  std::cout << "building cells in block...";
  build_cells_in_block();
  std::cout << "OK" << std::endl;

  build_support_regions();
}



void MultiScaleDataMSRSB::build_partitioning()
{
    auto & layer = active_layer();
    PureConnectionMap cell_connections;

    for (auto it = grid.begin_faces(); it != grid.end_faces(); ++it)
    {
      const auto & neighbors = it.neighbors();
      if (neighbors.size() == 2)  // not a boundary face
        cell_connections.insert_connection( neighbors[0], neighbors[1] );
    }

    layer.partitioning = multiscale::MetisInterface<hash_algorithms::empty>
        ::build_partitioning(cell_connections, layer.n_blocks, layer.n_cells);
}


void MultiScaleDataMSRSB::build_support_regions()
{
  std::cout << "finding block centroids...";
  find_centroids();
  std::cout << "OK" << std::endl;

  std::cout << "building block connections...";
  build_block_connections();
  std::cout << "OK" << std::endl;

  std::cout << "build block face data..." << std::flush;
  build_block_face_data();
  std::cout << "OK" << std::endl;

  std::cout << "build support regions..." << std::flush;
  active_layer().support_internal.resize(active_layer().n_blocks);
  for (std::size_t block = 0; block < active_layer().n_blocks; ++block)
  {
    std::cout << "\rbuild support regions... block "
              << (int)(100*block/active_layer().n_blocks)  << "% " <<std::flush;
    build_support_region(block);
  }

  std::cout << "OK" << std::endl;
}


void MultiScaleDataMSRSB::find_centroids()
{
  auto & layer = active_layer();
  layer.block_centroids.resize(layer.n_blocks);
  vector<size_t> n_cells_per_block(layer.n_blocks);

  for (auto cell = grid.begin_cells(); cell != grid.end_cells(); ++cell)
  {
    const size_t block = layer.partitioning[cell.index()];
    layer.block_centroids[block] += cell.center();
    n_cells_per_block[block]++;
  }

  for (std::size_t block=0; block<layer.n_blocks; ++block)
    layer.block_centroids[block] /= n_cells_per_block[block];

  // find closest cell
  layer.coarse_to_fine.resize(layer.n_blocks);
  const auto max = std::numeric_limits<double>::max();
  for (size_t block = 0; block < layer.n_blocks; ++block )
  {
    double min_dist = max;
    size_t closest = 0;
    for (const size_t cell : layer.cells_in_block[block])
    {
      const double current_dist = layer.block_centroids[block].distance(grid.get_center(cell));
      if ( current_dist <= min_dist )
      {
        closest = cell;
        min_dist = current_dist;
      }
    }
    layer.coarse_to_fine[block] = closest;
  }
}


void MultiScaleDataMSRSB::build_block_connections()
{
  auto & layer = active_layer();

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    if (neighbors.size() == 2)  // not a boundary face
    {
      const std::size_t i1 = layer.partitioning[neighbors[0]];
      const std::size_t i2 = layer.partitioning[neighbors[1]];;
      if (i1 != i2)
        if (!layer.block_internal_connections.connection_exists(i1, i2))
          layer.block_internal_connections.insert_connection(i1, i2);
    }
  }
}

void MultiScaleDataMSRSB::build_block_face_data()
{
  algorithms::UnionFindWrapper<size_t> face_disjoint = build_external_face_disjoint();

  size_t ghost_block = active_layer().n_blocks;
  // face group -> ghost block
  std::unordered_map<size_t, size_t> map_block_group;

  for ( const auto &it_face: face_disjoint.items() )
    if (map_block_group.find(face_disjoint.group( it_face.second )) == map_block_group.end())
      map_block_group.insert({ face_disjoint.group( it_face.second ), ghost_block++});

  find_block_face_centroids(face_disjoint, map_block_group);
  const auto map_vertex_blocks = build_map_vertex_blocks(face_disjoint, map_block_group);
  const auto map_block_edge_vertices = build_block_edges(map_vertex_blocks);
  find_block_edge_centroids(map_block_edge_vertices);
  debug_make_ghost_cell_names(face_disjoint, map_block_group);
}


void MultiScaleDataMSRSB::find_block_face_centroids(algorithms::UnionFindWrapper<size_t> & face_disjoint,
                                                    std::unordered_map<size_t, size_t>   & map_block_group)
{
  auto & layer = active_layer();
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    size_t block1, block2;
    block1 = layer.partitioning[neighbors[0]];

    if (neighbors.size() == 1)  // block + ghost block
    {
      const size_t face_group = face_disjoint.group(face.index());
      block2 = map_block_group[face_group];
    }
    else // if (neighbors.size() == 2)
    {
      block2 = layer.partitioning[neighbors[1]];
      if (block1 == block2) continue;
    }

    size_t conn_ind;
    if (layer.block_faces.connection_exists(block1, block2))
      conn_ind = layer.block_faces.connection_index(block1, block2);
    else
      conn_ind = layer.block_faces.insert_connection(block1, block2);

    auto & data = layer.block_faces.get_data(conn_ind);
    data.center += face.center();
    data.n_cell_faces++;
  }

  for (auto face = layer.block_faces.begin(); face != layer.block_faces.end(); ++face)
    face->center /= face->n_cell_faces;
}


bool MultiScaleDataMSRSB::share_edge(const mesh::const_face_iterator &face1,
                                     const mesh::const_face_iterator &face2)
{
  unordered_set<size_t> shared_verts;
  for (const auto &v1 : face1.vertex_indices())
    for (const auto &v2 : face2.vertex_indices())
      if (v1 == v2)
        shared_verts.insert(v1);

  return (shared_verts.size() > 1);
}


algorithms::UnionFindWrapper<size_t> MultiScaleDataMSRSB::build_external_face_disjoint()
{
  auto & layer = active_layer();

  // build map vertex - external face and fill disjoints
  std::unordered_map<size_t, vector<mesh::const_face_iterator>> map_vertex_face;
  algorithms::UnionFindWrapper<size_t> face_disjoint;

  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
    if (face.neighbors().size() == 1)
    {
      face_disjoint.insert(face.index());

      for (size_t & v : face.vertex_indices() )
      {
        auto it = map_vertex_face.find(v);
        if (it == map_vertex_face.end()) map_vertex_face.insert({v, {std::move( face )}});
        else it->second.push_back(std::move( face ));
      }
    }

  face_disjoint.finalize();

  //  create external face groups
  for (const auto & it : map_vertex_face)
  {
    const auto & faces = it.second;
    for (std::size_t i=0; i<faces.size(); ++i)
    {
      const auto & face1 = faces[i];
      const size_t cell1 = face1.neighbors()[0];
      for (std::size_t j=i+1; j<faces.size(); ++j)
      {
        const auto & face2 = faces[j];
        const size_t cell2 = face2.neighbors()[0];
        // if (face1 != face2) // they are created different

        /* criteria for uniting faces into groups to identify
         * ghost blocks:
         * (1) two boundary face do not belong to the  same cell
         * (2) their normals don't have a crazy angle between 'em
         * (3) two faces share an edge (two connected vertices) */
        if (cell1 != cell2)                                      // criterion 1
          if ( abs(dot(face1.normal(), face2.normal())) > 1e-4 ) // criterion 2
            if ( share_edge(face1, face2) )                      // criterion 3
            face_disjoint.merge(face1.index(), face2.index());
      }
    }
  }
  return face_disjoint;
}


std::unordered_map<std::size_t, std::vector<std::size_t>>
MultiScaleDataMSRSB::
build_map_vertex_blocks(algorithms::UnionFindWrapper<size_t> & face_disjoint,
                        std::unordered_map<size_t, size_t>   & map_block_group)
{
  /* collect vertices from faces that are on block-block interfaces */
  auto & layer = active_layer();
  std::unordered_map<std::size_t, std::vector<std::size_t>> vertex_blocks;
  for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
  {
    const auto & neighbors = face.neighbors();
    size_t block1, block2;
    block1 = layer.partitioning[neighbors[0]];

    if (neighbors.size() == 2)
    {
      block2 = layer.partitioning[neighbors[1]];
      if (block1 == block2) continue;
    }
    else if (neighbors.size() == 1)
    {
      const size_t face_group = face_disjoint.group(face.index());
      block2 = map_block_group[face_group];
    }

    for (const auto & vertex : face.vertex_indices())
    {
      auto it_vert = vertex_blocks.find(vertex);

      if (it_vert == vertex_blocks.end())
        vertex_blocks.insert({vertex, {block1, block2}});
      else
      {
        auto & blocks = it_vert->second;
        // std::find works here since a vertex probaby doesn't connect to
        // more than a few blocks
        if (std::find(blocks.begin(), blocks.end(), block1) == blocks.end())
          blocks.push_back(block1);
        if (std::find(blocks.begin(), blocks.end(), block2) == blocks.end())
          blocks.push_back(block2);
      }
    }
  }
  return vertex_blocks;
}


bool MultiScaleDataMSRSB::is_ghost_block(const size_t block) const
{
  if (block < active_layer().n_blocks)
    return false;
  else return true;
}


std::unordered_map<std::tuple<std::size_t,std::size_t,std::size_t>, std::vector<std::size_t>>
MultiScaleDataMSRSB::
build_block_edges(const std::unordered_map<std::size_t, std::vector<std::size_t>> & map_block_vertices)
{
  std::unordered_map<std::tuple<std::size_t,std::size_t,std::size_t>,
                     std::vector<std::size_t>> block_triplet_vertices;

  for (auto & it_vertex : map_block_vertices)
    if (it_vertex.second.size() > 2)  // on edge
    {
      // copy cause modifying (sort)
      auto blocks = it_vertex.second;
      // build block triplets with ascending ordering
      std::sort(blocks.begin(), blocks.end());

      for (std::size_t i = 0; i < blocks.size(); ++i)
        for (std::size_t j = i+1; j < blocks.size(); ++j)
          for (std::size_t k = j+1; k < blocks.size(); ++k)
          {
            const auto block_triplet = std::make_tuple(blocks[i], blocks[j], blocks[k]);
            auto it_blocks = block_triplet_vertices.find(block_triplet);
            if (it_blocks != block_triplet_vertices.end())
              it_blocks->second.push_back(it_vertex.first);
            else
              block_triplet_vertices.insert({std::move( block_triplet ), {it_vertex.first}});
          }
    }

  return block_triplet_vertices;
}


void MultiScaleDataMSRSB::find_block_edge_centroids(
      const std::unordered_map<std::tuple<std::size_t,std::size_t,std::size_t>, std::vector<std::size_t>>
      &map_block_edge_vertices)
{
  auto & layer = active_layer();

  for (const auto & block_edge : map_block_edge_vertices)
  {
    // find block edge center
    angem::Point<3,double> block_edge_center = {0, 0, 0};
    for (const auto & vertex : block_edge.second)
      block_edge_center += grid.vertex_coordinates(vertex);
    block_edge_center /= block_edge.second.size();

    // move this point to lie on the edge
    // make this point lie on edge
    block_edge_center = grid.vertex_coordinates(
        angem::find_closest_vertex(block_edge_center,
              /* all_vertices = */ grid.get_vertices(),
                    /* subset = */ block_edge.second));

    std::pair<std::size_t, std::size_t> face_blocks;
    std::size_t edging_block;
    for (std::size_t i = 0; i < 3; ++i)
    {
      std::size_t i1, i2;
      switch (i)
      {
        case 0:
          edging_block = std::get<0>(block_edge.first);
          i1 = std::get<1>(block_edge.first);
          i2 = std::get<2>(block_edge.first);
          break;
        case 1:
          edging_block = std::get<1>(block_edge.first);
          i1 = std::get<0>(block_edge.first);
          i2 = std::get<2>(block_edge.first);
          break;
        case 2:
          edging_block = std::get<2>(block_edge.first);
          i1 = std::get<0>(block_edge.first);
          i2 = std::get<1>(block_edge.first);
          break;
      }

      BlockFace * p_face;
      if (layer.block_faces.connection_exists(i1, i2))
        p_face = &(layer.block_faces.get_data(i1, i2));
      else  // connection only by edge
      {
        if (is_ghost_block(i1) or is_ghost_block(i2)) continue;

        const std::size_t con = layer.block_faces.insert_connection(i1, i2);
        p_face = &(layer.block_faces.get_data(con));
        p_face->center = block_edge_center;
      }

      auto & face = *p_face;
      if (std::find(face.edging_blocks.begin(), face.edging_blocks.end(),
                    edging_block) == face.edging_blocks.end())
      {
        face.edging_blocks.push_back(edging_block);
        face.edge_centers.push_back(block_edge_center);
      }
    }  // end triplet combination loop
  }    // end block loop
}


bool have_better_fit(const std::size_t                                 block,
                     const std::size_t                                 n1,
                     const std::size_t                                 n2,
                     const PureConnectionMap                         & block_connections,
                     const hash_algorithms::ConnectionMap<BlockFace> & block_faces)
{
  for (const size_t nn1 : block_faces.get_neighbors(n1))
  {
    if (block_connections.connection_exists(nn1, n1) &&
        block_connections.connection_exists(nn1, n2))
    {
      // std::cout << "better = " << nn1 <<" ";
      return true;
    }

    if (block_connections.connection_exists(nn1, n1) &&
        block_faces.connection_exists(nn1, n2))
      return have_better_fit(block, nn1, n2, block_connections, block_faces);
  }
  return false;
}


void MultiScaleDataMSRSB::build_support_region(const std::size_t block)
{
  /* This algorithm is geometric and by far does not work in all cases.
   * It at least requires some levels of convexity of the coarse blocks
   * though a rigorous convexity is not required.
   * The basic idea is to select a surface comprised of triangles
   * that limits the support region.
   * All the cells and vertices limited by such planes are within
   * the support region.
   */
  auto & layer = active_layer();
  mesh::SurfaceMesh<double> bounding_surface;

  // exit(0);
  // std::cout << "\nblock = " << block << std::endl;
  // for (auto it = layer.block_faces.begin(); it != layer.block_faces.end(); ++it)
  //   std::cout << it.elements().first <<" " << it.elements().second << std::endl;

  //  select non-ghost neighobors of the current block
  for (const size_t & neighbor1 : layer.block_faces.get_neighbors(block) )
    if (!is_ghost_block(neighbor1))
    {
      // std::cout << "\tneighbor1 = " << neighbor1 << std::endl;
      // select neighbors of neighbor1 (including ghosts)
      for (const size_t & neighbor2 : layer.block_faces.get_neighbors(neighbor1))
      {
        // std::cout << "\t\tneighbor2 = " << debug_ghost_cell_names[ neighbor2 ] <<"...";

        // criterion 1: should be a neighbor of the block
        if (!layer.block_faces.connection_exists(neighbor2, block))
        {
          // std::cout << "no conn with " << block << std::endl;
          // std::cout << "-" << std::endl;
          continue;
        }

        // criterion 2: the face between neighbor2 and block should have a
        // common edge with neighbor1
        const auto & b_n2_edging = layer.block_faces.get_data(neighbor2, block).edging_blocks;
        // Note: we use std::find since the number of edging blocks is assumed <O(10)
        if (find( b_n2_edging.begin(), b_n2_edging.end(), neighbor1 ) == b_n2_edging.end())
        {
          // std::cout << neighbor1 << " not in " << block <<"-"<<neighbor2 <<" edging" << std::endl;
          continue;
        }

        // critetion 3: treatment of multiple blocks having a common coarse edge
        // in this case we select the neighbors that actually have a common face
        // this is done recursively with have_better_fit
        // NOTE: of course those are not ghost blocks
        if (!layer.block_internal_connections.connection_exists(neighbor1, neighbor2)) // no conn by face
          if (!is_ghost_block(neighbor1) and !is_ghost_block(neighbor2))
            if (have_better_fit(block, neighbor1, neighbor2,
                                layer.block_internal_connections,
                                layer.block_faces))
            {
              // std::cout << "have better fit" << std::endl;
              continue;
            }

        // select third neighbor to determine the direction of the plane that
        // limits the support region
        const auto & face_n1_n2 = layer.block_faces.get_data(neighbor1, neighbor2);
        for (const size_t & neighbor3 : face_n1_n2.edging_blocks)
        {
          // std::cout << "\n\t\t\tneighbor3 = " << debug_ghost_cell_names[ neighbor3 ] <<" ";
          // Obviously this should'n be an already selected block
          if (neighbor3 == neighbor1 or neighbor3 == neighbor2 or neighbor3 == block)
          {
            // std::cout << "repeated block" << std::endl;
            continue;
          }

          // Criterion 1 (yes, again)
          if (!layer.block_faces.connection_exists(neighbor3, block))
          {
            // std::cout << "not connected to " << block << " face!"<< std::endl;
            continue;
          }
          // and again
          if (!layer.block_faces.connection_exists(neighbor1, neighbor3))
          {
            // std::cout << "not connected to "<< neighbor1 << std::endl;
            continue;
          }

          // Criterion 2 (same)
          const auto & b_n3_edging = layer.block_faces.get_data(neighbor3, block).edging_blocks;
          if (find( b_n3_edging.begin(), b_n3_edging.end(), neighbor2 ) == b_n3_edging.end())
          {
            // std::cout << neighbor2 <<" not edging in "<< block<<"-"<<neighbor3<< std::endl;
            continue;
          }
          // Criterion 2 (same) for another pair of blocks
          if (find( b_n3_edging.begin(), b_n3_edging.end(), neighbor1 ) == b_n3_edging.end())
          {
            // std::cout << neighbor1 <<" not edging in "<< block<<"-"<<neighbor3<< std::endl;
            continue;
          }

          // Criterion 3 (same)
          if (!layer.block_internal_connections.connection_exists(neighbor1, neighbor3)) // no conn by face
            if (!is_ghost_block(neighbor1) and !is_ghost_block(neighbor3))
              if (have_better_fit(block, neighbor1, neighbor3,
                        layer.block_internal_connections, layer.block_faces))
              {
                // std::cout << "better fit" << std::endl;
                continue;
              }
          // and again
          if (!layer.block_internal_connections.connection_exists(neighbor2, neighbor3)) // no conn by face
            if (!is_ghost_block(neighbor2) and !is_ghost_block(neighbor3))
              if (have_better_fit(block, neighbor2, neighbor3,
                        layer.block_internal_connections, layer.block_faces))
              {
                // std::cout << "better fit" << std::endl;
                continue;
              }

          // std::cout << "OK" << std::endl;

          // Now it's time to construct an element of the bounding surface
          //  construct a triangle: a part of the bounding surface
          // const std::size_t index_of_n3 = find_index(face_n1_n2.edging_blocks, neighbor3);
          const std::size_t index_of_n3 = find(face_n1_n2.edging_blocks.begin(),
                                               face_n1_n2.edging_blocks.end(),
                                               neighbor3) -
                                          face_n1_n2.edging_blocks.begin();

          const std::vector<Point> bounding_triangle_vertices =
              {layer.block_centroids[neighbor1], face_n1_n2.center,
               face_n1_n2.edge_centers[index_of_n3]};

          angem::Polygon<double> bounding_triangle(bounding_triangle_vertices);
          // store triangle for later
          bounding_surface.insert(bounding_triangle);
          build_support_region_boundary(block, neighbor1, bounding_triangle);
        }
      }
    }

  // debug output bounding surface vtk
  // const std::string fname = "debug_output/support_surface-" + std::to_string(block) + ".vtk";
  // IO::VTKWriter::write_surface_geometry(bounding_surface.get_vertices(),
  //                                       bounding_surface.get_polygons(), fname);

  // Next we find the internal points of the support region
  build_support_internal_cells(block, bounding_surface);
}


void MultiScaleDataMSRSB::build_support_region_boundary(const std::size_t block,
                                                        const std::size_t neighbor,
                                                        const angem::Shape<double> & bounding_shape)
{
  auto & layer = active_layer();
  if (layer.support_boundary.empty())
    layer.support_boundary.resize(layer.n_blocks);

  // fast collision checking algorithm
  angem::CollisionGJK<double> collision;
  for (const std::size_t cell : layer.cells_in_block[neighbor])
  {
    const auto p_cell_polyhedra = grid.get_polyhedron(cell);
    if (collision.check(*p_cell_polyhedra, bounding_shape))
      layer.support_boundary[block].insert(cell);
  }
}


Point MultiScaleDataMSRSB::find_point_outside_support_region(const std::size_t block)
{
  /* return the center of the block that is does not have a face with the current block */
  const auto & layer = active_layer();
  const auto & block_neighbors = layer.block_internal_connections.get_neighbors(block);

  for (std::size_t i = 0 ; i < layer.n_blocks; ++i)
    if (block != i)
      if(find(block_neighbors.begin(), block_neighbors.end(), i) == block_neighbors.end())
        return layer.block_centroids[i];

  //  just to shut the compiler up
  return {1e7, 1e7, 1e7};
}


void MultiScaleDataMSRSB::build_support_internal_cells(const std::size_t block,
                                                       const mesh::SurfaceMesh<double>& bounding_surface)
{
  auto & layer = active_layer();

  for (const std::size_t cell: layer.cells_in_block[block])
    layer.support_internal[block].insert(cell);

  // loop through block neighbours and determine which cells
  // belong to the support region
  // ray casting algorithm:
  // if a ray starting from the outer point and going into the
  // cell center is intersecting the bounding surface an odd number
  // of times, then the cell center is inside the support region
  const Point outer_point = find_point_outside_support_region(block);
  const auto & block_neighbors = layer.block_internal_connections.get_neighbors(block);
  // std::cout << "outer_point = " << outer_point << std::endl;

  for (const std::size_t neighbor : block_neighbors)
    for (const std::size_t cell : layer.cells_in_block[neighbor])
    {
      // if (block == 1) std::cout << cell << ": " << grid.get_center(cell) << " ";
      // ray casting
      if ( angem::point_inside_surface(grid.get_center(cell), outer_point,
                                       bounding_surface.get_vertices(),
                                       bounding_surface.get_polygons()) )
      {
        // if (block == 1) std::cout << "OK" << std::endl;
        // if not among the boundary cells
        if (layer.support_boundary[block].find(cell) == layer.support_boundary[block].end())
          layer.support_internal[block].insert(cell);
      }
      // else if (block == 1) std::cout << "NO" << std::endl;
    }
}


void MultiScaleDataMSRSB::fill_output_model(MultiScaleOutputData & model,
                                            const int layer_index) const
{
  const auto & layer = layers[layer_index];

  model.cell_data = true;
  model.n_coarse = layer.n_blocks;

  // partitioning
  model.partitioning.resize(layer.partitioning.size());
  std::copy(layer.partitioning.begin(), layer.partitioning.end(), model.partitioning.begin());

  //  centroids
  model.centroids = layer.coarse_to_fine;

  // support boundary
  model.support_boundary = layer.support_boundary;

  //  support internal
  model.support_internal = layer.support_internal;
}


void MultiScaleDataMSRSB::build_cells_in_block()
{
  auto & layer = active_layer();
  layer.cells_in_block.resize(layer.n_blocks);
  for (std::size_t cell=0; cell<layer.n_cells; ++cell)
    layer.cells_in_block[layer.partitioning[cell]].push_back(cell);
}

std::string debug_direction_to_string(angem::Point<3,double> p)
{
  angem::Point<3,double> pn = p; pn.normalize();
  // std::cout << "pn = " << pn << std::endl;
  double tol = 1e-4;
  repeat:{
    if ((pn - angem::Point<3,double>(-1, 0, 0)).norm() < tol)
      return "left";
    if ((pn - angem::Point<3,double>(1, 0, 0)).norm() < tol)
      return "right";
    if ((pn - angem::Point<3,double>(0, -1, 0)).norm() < tol)
      return "front";
    if ((pn - angem::Point<3,double>(0, 1, 0)).norm() < tol)
      return "back";

    if ((pn - angem::Point<3,double>(0, 0, 1)).norm() < tol)
      return "top";
    if ((pn - angem::Point<3,double>(0, 0, -1)).norm() < tol)
      return "bottom";
    tol *= 2;
    goto repeat;
  }

  return "unknown";
}

void MultiScaleDataMSRSB::
debug_make_ghost_cell_names(const algorithms::UnionFindWrapper<size_t> & face_disjoint,
                            std::unordered_map<size_t, size_t>   & map_block_group)
{
  const auto & layer = active_layer();
  // std::vector<Point> avg_face_group_direction(face_disjoint.n_groups());
  // std::vector<size_t> n_group_items(face_disjoint.n_groups());

  // Point global_center;
  // size_t n_items;
  // std::cout << "starting loop" << std::endl;

  for (auto & group : map_block_group)
  {
    const auto ghost = group.second;
    // std::cout << "ghost = " << ghost << std::endl;
    const auto block = layer.block_faces.get_neighbors(ghost)[0];
    // std::cout << "block = " << block << std::endl;
    if (block < layer.n_blocks)
    {
      Point face_direction;
      size_t n_faces = 0;
      // find ghost face center
      for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
        if (face.neighbors().size() == 1)
          if (layer.partitioning[face.neighbors()[0]] == block)
            if ( face_disjoint.group(face.index()) == group.first)
          {
            face_direction += face.center();
            n_faces++;
          }

      face_direction /= n_faces;
      const std::string block_name = debug_direction_to_string(face_direction -
                                                               layer.block_centroids[block]);
      // std::cout << "adding ghost cell  " << block_name  << std::endl;
      debug_ghost_cell_names[ghost] = block_name;

    }
  }

  for (size_t i = 0; i < layer.n_blocks; i++)
    debug_ghost_cell_names[i] = std::to_string(i);

  // // print which blocks connected to which ghosts
  // for (size_t i = 0; i < layer.n_blocks; i++)
  // {
  //   std::cout << "block " << i <<": ";
  //   for (const auto neighbor : layer.block_faces.get_neighbors(i))
  //     if (neighbor >= layer.n_blocks)
  //     {
  //       std::cout << debug_ghost_cell_names[neighbor] << "\t";
  //     }
  //   std::cout << std::endl;
  // }
}

}  // end namespace
