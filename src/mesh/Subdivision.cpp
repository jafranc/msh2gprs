#include "Subdivision.hpp"

namespace mesh {

using Point = angem::Point<3,double>;

Subdivision::Subdivision(const mesh::Cell & cell, mesh::Mesh & triangulation, const size_t order)
    : _parent_cell(cell), _grid(triangulation), _order(order)
{
  if (!_grid.empty())
    throw std::invalid_argument("Grid should be empty");

  // copy cell from its master grid to a new grid
  // this cell will be the parent of all cells in the grid
  create_master_cell_();

  perform_subdivision_r0_(*_grid.begin_active_cells());
  for (size_t i=0; i<order; ++i)
  {
    std::vector<size_t> active_cells(_grid.n_active_cells());
    size_t j = 0;
    for (auto cell = _grid.begin_active_cells(); cell != _grid.end_active_cells(); ++cell, ++j)
      active_cells[j] = cell->index();
    for (const size_t icell : active_cells)
      perform_subdivision_tetra_(_grid.cell(icell));
  }
}

void Subdivision::create_master_cell_()
{
  _grid.vertices() = _parent_cell.vertex_coordinates();
  std::map<size_t,size_t> old_to_new;
  size_t i = 0;
  for (const size_t v : _parent_cell.vertices())
    old_to_new[v] = i++;
  const auto faces = _parent_cell.faces();
  std::vector<FaceTmpData> new_faces(faces.size());
  i = 0;
  for (const auto face : faces)
  {
    const std::vector<std::size_t> & vertices = face->vertices();
    auto & f = new_faces[i];
    f.vertices.reserve(vertices.size());
    for (size_t v=0; v<vertices.size(); ++v)
      f.vertices.push_back( old_to_new[vertices[v]] );
    f.vtk_id = face->vtk_id();
    f.marker = i + 1;  // to comply with old face markering by gmsh
    i++;
  }

  std::vector<size_t> use_faces(new_faces.size());
  std::iota( use_faces.begin(), use_faces.end(), 0 );
  std::vector<size_t> ivertices(_grid.n_vertices());
  std::iota( ivertices.begin(), ivertices.end(), 0 );
  // std::cout << std::flush << std::endl;
  _grid.insert_cell_(ivertices, use_faces, new_faces,
                     _parent_cell.vtk_id(), _parent_cell.marker());
}

void Subdivision::perform_subdivision_r0_(mesh::Cell & cell)
{
  const size_t parent_cell_index = cell.index();
  const size_t cell_center_index = _grid.n_vertices();
  assert( cell.is_active() );
  _grid.vertices().push_back( cell.center() );
  // save face indices since inserting new cells invalidates pointsers
  std::vector<size_t> face_indices;
  for( auto face :cell.faces() )
    face_indices.push_back( face->index() );

  for (const size_t iface : face_indices)
  {
    const mesh::Face * face = &_grid.face(iface);
    const auto c = face->center();
    size_t face_center_index = _created_vertices.find(c);
    if ( face_center_index == _created_vertices.size() )
    {
      _created_vertices.insert(c);
      _created_vertex_indices.push_back(_grid.insert_vertex(c));
      face_center_index = _created_vertex_indices.back();
    }

    auto build_trgl_face = [](const std::vector<size_t> vertices, const size_t parent,
                              const int marker, FaceTmpData & f) {
                             f.vertices = vertices;
                             f.vtk_id = angem::VTK_ID::TriangleID;
                             f.parent = parent;
                             f.marker = marker;
                           };

    // refine face and build tetras
    const auto face_vertices = face->vertices();
    for (size_t i=0; i<face_vertices.size(); ++i)
    {
      size_t i1 = face_vertices[i], i2 = face_vertices[i+1];
      if (i == face_vertices.size() - 1) i2 = face_vertices[0];
      std::vector<FaceTmpData> tetra_faces(4);
      // // base of the tetra resides on the parent face
      build_trgl_face({i1, i2, face_center_index}, face->index(), face->marker(), tetra_faces[0]);
      // add three more faces to build the tetrahedron
      build_trgl_face({i1, i2, cell_center_index}, constants::invalid_index,
                      constants::default_face_marker, tetra_faces[1]);
      build_trgl_face({i1, face_center_index, cell_center_index}, constants::invalid_index,
                      constants::default_face_marker, tetra_faces[2]);
      build_trgl_face({i2, face_center_index, cell_center_index}, constants::invalid_index,
                      constants::default_face_marker, tetra_faces[3]);
      std::vector<size_t> take_faces(4);
      std::iota(take_faces.begin(), take_faces.end(), 0);
      const size_t child_cell_index =
          _grid.insert_cell_({i1, i2, face_center_index, cell_center_index},
                             take_faces, tetra_faces, angem::TetrahedronID,
                             cell.marker());
      // const size_t child_cell_index =
      // _grid.insert_cell( {i1, i2, face_center_index, cell_center_index},
      //                    angem::TetrahedronID, cell.marker());
      _grid.m_cells[parent_cell_index].m_children.push_back(child_cell_index);
      _grid.m_cells[child_cell_index].m_parent = cell.index();
    }
  }
  _grid.m_n_cells_with_hanging_nodes++;
}

void Subdivision::perform_subdivision_tetra_(Cell & cell)
{
  /*
   * Bey, Jürgen. "Tetrahedral grid refinement." Computing 55.4 (1995): 355-378.
   *
   *                         Illustration
   *             2                                     2
   *           ,/|`\                                 ,/|`\
   *         ,/  |  `\                             ,/  |  `\
   *       ,/    '.   `\                         ,6    '.   `5
   *     ,/       |     `\                     ,/       8     `\
   *   ,/         |       `\                 ,/         |       `\
   *  0-----------'.--------1               0--------4--'.--------1
   *   `\.         |      ,/                 `\.         |      ,/
   *      `\.      |    ,/                      `\.      |    ,9
   *         `\.   '. ,/                           `7.   '. ,/
   *            `\. |/                                `\. |/
   *               `3                                    `3
   */

  assert( cell.is_active() );
  const auto vertex_coord = cell.vertex_coordinates();
  std::vector<Point> new_vertex_coord = {0.5 * (vertex_coord[0]  + vertex_coord[1]),  // 4
                                         0.5 * (vertex_coord[2]  + vertex_coord[1]),  // 5
                                         0.5 * (vertex_coord[0]  + vertex_coord[2]),  // 6
                                         0.5 * (vertex_coord[0]  + vertex_coord[3]),  // 7
                                         0.5 * (vertex_coord[2]  + vertex_coord[3]),  // 8
                                         0.5 * (vertex_coord[1]  + vertex_coord[3])}; // 9
  std::vector<size_t> vertices = cell.vertices();
  for (auto & p : new_vertex_coord)
  {
    const size_t idx = _created_vertices.find(p);
    if ( idx == _created_vertices.size() )  // new vertex
    {
      _created_vertices.insert(p);
      _created_vertex_indices.push_back(_grid.insert_vertex(p));
      vertices.push_back(_created_vertex_indices.back());
    }
    else  // vertex exists
      vertices.push_back(_created_vertex_indices[idx]);
  }

  // Tetra 1
  const size_t parent_cell_index = cell.index();
  insert_tetra_({0, 4, 6, 7}, vertices,
                {{0, 4, 6}, {0, 4, 7},
                 {0, 6, 7}, {4, 6, 7}}, parent_cell_index);
  // Tetra 2
  insert_tetra_({2, 5, 6, 8}, vertices,
                {{2, 6, 8}, {2, 5, 8},
                 {2, 5, 6}, {5, 6, 8}}, parent_cell_index);
  // Tetra 3
  insert_tetra_({4, 6, 7, 8}, vertices,
                {{4, 6, 7}, {4, 7, 8},
                 {4, 6, 8}, {6, 7, 8}}, parent_cell_index);
  // Tetra 4
  insert_tetra_({4, 5, 6, 8}, vertices,
                {{4, 5, 6}, {4, 5, 8},
                 {4, 6, 8}, {5, 6, 8}}, parent_cell_index);
  // tetra 5
  insert_tetra_({3, 7, 8, 9}, vertices,
                {{3, 7, 8}, {3, 8, 9},
                 {3, 7, 9}, {7, 8, 9}}, parent_cell_index);
  // Tetra 6
  insert_tetra_({4, 7, 8, 9}, vertices,
                {{4, 7, 8}, {4, 7, 9},
                 {4, 8, 9}, {7, 8, 9}}, parent_cell_index);
  // Tetra 7
  insert_tetra_({4, 5, 8, 9}, vertices,
                {{4, 5, 8}, {4, 5, 9},
                 {4, 8, 9}, {5, 8, 9}}, parent_cell_index);
  // Tetra 8
  insert_tetra_({1, 4, 5, 9}, vertices,
                {{1, 4, 5}, {1, 4, 9},
                 {1, 5, 9}, {4, 5, 9}}, parent_cell_index);
  _grid.m_n_cells_with_hanging_nodes++;
}

void Subdivision::insert_tetra_(const std::vector<size_t> & local_vertex_indices,
                                const std::vector<size_t> & global_vertex_indices,
                                const std::vector<std::vector<size_t>> & faces,
                                const size_t parent_cell_index)
{
  std::vector<size_t> cell_vertices(local_vertex_indices.size());
  for (size_t v=0; v<local_vertex_indices.size(); ++v)
    cell_vertices[v] = global_vertex_indices[ local_vertex_indices[v] ];

  std::vector<FaceTmpData> cell_faces(4);
  for (size_t f=0; f<faces.size(); ++f)
  {
    auto & face = cell_faces[f];
    face.vertices.resize(3);
    face.vtk_id = angem::TriangleID;
    std::transform(faces[f].begin(), faces[f].end(), face.vertices.begin(),
                   [global_vertex_indices](size_t v) {return global_vertex_indices[v];});
  }
  std::vector<size_t> take_faces(4);
  std::iota(take_faces.begin(), take_faces.end(), 0);
  const size_t new_cell_index =
      _grid.insert_cell_(cell_vertices, take_faces, cell_faces, angem::TetrahedronID,
                         _parent_cell.marker());

  _grid.m_cells[parent_cell_index].m_children.push_back(new_cell_index);
  _grid.m_cells[new_cell_index].m_parent = parent_cell_index;
}


}  // end namespace mesh
