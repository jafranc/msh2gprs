#include "VTKReader.hpp"
#include <set>  // provides set

namespace mesh {

namespace io {

VTKReader::VTKReader(const std::string file_name, mesh::Mesh & grid)
    : _grid(grid)
{
  read_file_(file_name);
}

void VTKReader::read_file_(const std::string & file_name)
{
  std::fstream in;
  in.open(file_name.c_str(), std::fstream::in);
  if (!in) throw std::out_of_range(file_name + " does not exist");

  read_header_(in);
  read_vertices_(in);
  read_cells_(in);
  read_cell_types_(in);
  create_grid_();
  read_data_arrays_(in);
  in.close();
}

void VTKReader::read_header_(std::fstream & in) const
{
  std::string line;
  std::getline(in, line);  // vtk DataFile Version 3.0
  std::getline(in, line);  // 3D Grid
  std::getline(in, line);  // ASCII
  if (line != "ASCII")
    throw std::invalid_argument("File format not supported");
  std::getline(in, line);  // DATASET UNSTRUCTURED_GRID
}

void VTKReader::read_vertices_(std::fstream & in)
{
  std::string s;
  in >> s;
  if (s != "POINTS")
    throw std::invalid_argument("Vertices should start with POINTS");

  size_t n_vertices;
  in >> n_vertices;
  if (n_vertices == 0)
    throw std::invalid_argument("POINTS should have more than 0 vertices");
  std::cout << "n_vertices = " << n_vertices << std::endl;

  in >> s;
  if (s != "float")
    throw std::invalid_argument("Unknown POINTS format (should be FLOAT)");

  auto & vertex_coord = _grid.vertices();
  vertex_coord.resize(n_vertices);
  for (size_t i=0; i<n_vertices; ++i)
    for (size_t j=0; j<3; ++j)
      in >> vertex_coord[i][j];
}

void VTKReader::read_cells_(std::fstream & in)
{
  std::string s;
  in >> s;
  if ( s != "CELLS" )
    throw std::invalid_argument("Entry should be CELLS");
  size_t n_cells;
  in >> n_cells;
  std::cout << "n_cells = " << n_cells << std::endl;

  size_t n_cell_entries;
  in >> n_cell_entries;
  _cell_entries.resize(n_cell_entries);

  for (size_t i=0; i<n_cell_entries; ++i)
    in >> _cell_entries[i];
}

void VTKReader::read_cell_types_(std::fstream & in)
{
  std::string s;
  in >> s;
  if (s != "CELL_TYPES")
  {
    std::cout << s << std::endl;
    throw std::invalid_argument("Entry should be CELL_TYPES");
  }

  size_t n_cell_type_entries;
  in >> n_cell_type_entries;
  _vtk_ids.resize(n_cell_type_entries);
  for (size_t i=0; i<n_cell_type_entries; ++i)
    in >> _vtk_ids[i];
}

void VTKReader::create_grid_()
{
  const size_t n_cells = _vtk_ids.size();
  size_t idx = 0;

  //check if 2D mesh upgrade it to 3D by single layer extrusion
  bool is2D = true;
  std::for_each(_vtk_ids.cbegin(),_vtk_ids.cend(),[&is2D](const int ids){ is2D &=(ids<10);});
  if (is2D)
  {
    const size_t off_vert = _grid.vertices().size();
    _grid.vertices().resize( _grid.vertices().size() + off_vert);
    for( size_t icell = 0; icell < n_cells; ++icell )
    {
      const int id = _vtk_ids[icell];
      _vtk_ids[icell] = VTK_ID::GeneralPolyhedronID;
      extrude_to_gen_polyhedron(id, idx, off_vert);
    }
  }
  else
  {
    for( size_t icell = 0; icell < n_cells; ++icell )
    {
      const int id = _vtk_ids[icell];
      if( id == angem::VTK_ID::GeneralPolyhedronID )
        create_general_polyhedron_cell_( idx );
      else
        create_regular_polyhedron_cell_( id, idx );
    }
  }
}

void VTKReader::create_general_polyhedron_cell_(size_t & idx)
{
  const size_t n_cell_entries = _cell_entries[idx++];
  const size_t n_faces = _cell_entries[idx++];
  std::vector<mesh::FaceTmpData> faces(n_faces);
  for (size_t iface=0; iface<n_faces; ++iface)
  {
    const size_t nv = _cell_entries[idx++];
    auto & face = faces[iface];
    face.vertices.resize(nv);
    for (size_t i=0; i<nv; ++i)
      face.vertices[i] = _cell_entries[idx++];
  }

  std::vector<size_t> take_faces(n_faces);
  std::iota(take_faces.begin(), take_faces.end(), 0);
  std::set<size_t> cell_vertices_set;
  for (auto & f : faces)
    for (const size_t v : f.vertices)
      cell_vertices_set.insert(v);
  std::vector<size_t> cell_vertices(cell_vertices_set.begin(), cell_vertices_set.end());

  _grid.insert_cell_( cell_vertices, take_faces, faces , angem::VTK_ID::GeneralPolyhedronID,
                      mesh::constants::default_cell_marker);
}

void VTKReader::extrude_to_gen_polyhedron(const int id, size_t& idx, size_t off_vert )
{
  const size_t n_cell_entries = _cell_entries[idx++];
  const size_t n_faces = n_cell_entries + 2;


  std::vector<mesh::FaceTmpData> faces(n_faces);
  {
    const size_t nv = n_cell_entries;
    auto& face0 = faces.front();
    auto& facen = faces.back();
    face0.vertices.resize(nv);
    face0.vtk_id = id;
    facen.vertices.resize(nv);
    facen.vtk_id = id;


    for (size_t i=0; i<nv; ++i)
    {
      face0.vertices[i] = _cell_entries[idx];
      facen.vertices[i] = _cell_entries[idx++] + off_vert;

      //adding the proper vertices
      auto ncoords = _grid.vertices()[face0.vertices[i]];
      ncoords[2] += 1.0;
      _grid.m_vertices[facen.vertices[i]] = ncoords;
    }

    size_t iface = 1;
    const size_t base_sz = face0.vertices.size();
    for (size_t in=0; in < base_sz; ++in)
    {
      faces[iface].vertices.resize(4);
      faces[iface].vtk_id = 9;
      //using peridoc index to close the cycles
      faces[iface++].vertices = { face0.vertices[in], face0.vertices[((in+1)%base_sz)],
                                 facen.vertices[((in+1)%base_sz)], facen.vertices[in] };
    }
  }

  std::vector<size_t> take_faces(n_faces);
  std::iota(take_faces.begin(), take_faces.end(), 0);
  std::set<size_t> cell_vertices_set;
  for (auto & f : faces)
    for (const size_t v : f.vertices)
      cell_vertices_set.insert(v);
  std::vector<size_t> cell_vertices(cell_vertices_set.begin(), cell_vertices_set.end());


 _grid.insert_cell_( cell_vertices, take_faces, faces , angem::VTK_ID::GeneralPolyhedronID,
                     0 );
}


void VTKReader::create_regular_polyhedron_cell_(const int id, size_t & idx)
{
  const size_t n_cell_entries = _cell_entries[idx++];
  std::vector<size_t> cell_vertices(n_cell_entries);
  for (size_t i=0; i<n_cell_entries; ++i)
    cell_vertices[i] = _cell_entries[idx++];
  _grid.insert_cell(cell_vertices, id);
}

enum ArrayType
{
  undefined, cell_data, point_data
};

void VTKReader::read_data_arrays_(std::fstream & in)
{
  ArrayType current_section = undefined;
  while (!in.eof())
  {
    std::string array_type;
    in >> array_type;
    if (array_type == "CELL_DATA")
    {
      std::string tmp;
      in >> tmp;
      if (std::atoi(tmp.c_str()) != _grid.n_active_cells())
        throw std::runtime_error("Wrong size in CELL_DATA");
      in >> tmp;  // FIELD
      in >> tmp;  // FieldData
      size_t n_arrays;
      in >> n_arrays;
      for (size_t i = 0; i < n_arrays; ++i)
      {
        _cell_data.emplace_back();
        _cell_data_names.emplace_back();
        read_array_(in, _cell_data.back(), _cell_data_names.back(),_grid.n_active_cells());
      }
    }
    else if (array_type == "POINT_DATA")
    {
      std::string tmp;
      in >> tmp;
      if (std::atoi(tmp.c_str()) != _grid.n_vertices())
        throw std::runtime_error("Wrong size in CELL_DATA");
      in >> tmp;  // FIELD
      in >> tmp;  // FieldData
      _point_data_names.emplace_back();
      _point_data.emplace_back();
      read_array_(in, _point_data.back(), _point_data_names.back(), _grid.n_vertices());
    }
    else
    {
      if (array_type.empty())
        break;
    }
  }
}

void VTKReader::read_array_(std::fstream & in, std::vector<double> & data,
                            std::string & name, const size_t array_size) const
{
  if (data.size() != array_size)
    data.resize(array_size);
  in >> name;
  size_t n_comp; in >> n_comp;
  assert( n_comp == 1 );
  size_t n_entries; in >> n_entries;
  if (array_size != n_entries)
    throw std::runtime_error("invalid vtk file");
  std::string datatype; in >> datatype;
  size_t i = 0;
  while (!in.eof())
  {
    in >> data[i++];
    if (i == array_size)
      break;
  }
  if (i != array_size)
    throw std::runtime_error("invalid vtk file");
}

}  // end namespace io

}  // end namespace mesh
