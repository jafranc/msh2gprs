#pragma once

#include "angem/VTK_ID.hpp"
#include <vector>

namespace gprs_data {

using angem::VTK_ID;

enum GmshElementType
{
  node,
  edge,
  face,
  cell,
  invalid_element
};

// map gmsh element id to vtk id
extern std::vector<int> msh_id_to_vtk_id;
// entity of gmshs elements (edge, cell, node)
extern std::vector<GmshElementType> gmsh_element_types;
// number of vertices for each gmsh element
extern std::vector<std::size_t> gmsh_element_nvertices;

//2d to 3d
extern std::vector<int> msh_2d_3d_vtk;
extern std::vector<GmshElementType> msh_2d_3d;
extern std::vector<std::size_t> gmsh_element_nvert_2d_3d;


}  // end namespace gprs_data
