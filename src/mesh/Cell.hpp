#pragma once

#include "angem/Point.hpp"
#include "angem/Polyhedron.hpp"
#include "Face.hpp"
#include "constants.hpp"
#include <vector>
#include <functional>  // std::reference_wrapper

namespace mesh
{
using Point = angem::Point<3,double>;
using Polyhedron = angem::Polyhedron<double>;

class Cell
{
 public:
  // constructor
  Cell(const std::size_t                cell_index,
       const std::vector<std::size_t> & vertices,
       const std::vector<std::size_t> & faces,
       std::vector<Point>             & grid_vertices,
       std::vector<Cell>              & grid_cells,
       std::vector<Face>              & grid_faces,
       const int                        vtk_id,
       const int                        marker,
       const std::size_t                parent = constants::invalid_index);
  // assignment operator
  Cell & operator=(const Cell & other);
  // comparison operator
  inline bool operator==(const Cell & other) const { return index() == other.index(); }
  // comparison operator
  inline bool operator!=(const Cell & other) const { return !(*this == other); }
  // ---------------------- ACCESS OPERATORS ------------------------------- //
  // get cell index
  inline std::size_t index() const { return m_index; }
  // get cell marker
  inline int marker() const { return m_marker; }
  // get vtk id
  inline int vtk_id() const { return m_vtk_id; }
  // get vector of vertex indices
  inline std::vector<std::size_t> & vertices() { return m_vertices; }
  // get const vector of vertex indices
  inline const std::vector<std::size_t> & vertices() const { return m_vertices; }
  // get vector of neighbors
  std::vector<Cell*> neighbors();
  // get vector of neighbors
  std::vector<const Cell*> neighbors() const;
  // get vector of cell faces
  std::vector<Face*> faces();
  // get vector of cell faces
  std::vector<const Face*> faces() const;
  // convenience methods
  // get the number of cell vertices
  inline std::size_t n_vertices() const { return m_vertices.size(); }
  // get center of mass
  Point center() const;
  // get cell volume
  double volume() const;
  // get a polyhedron that represents a cell
  std::unique_ptr<Polyhedron> polyhedron() const;
  // true if cell hace a vertex, false otherwise
  bool has_vertex(const std::size_t vertex_index) const;
  // returns true if has no children; else returns false
  inline bool is_active() const {return m_children.empty();}
  // returns the parent cell. If cell has not parents, returns itself.
  inline const Cell & parent() const { return (*pm_grid_cells)[m_parent]; }
  // returns the parent cell. If cell has not parents, returns itself.
  inline Cell & parent() { return (*pm_grid_cells)[m_parent]; }
  // returns the parent of parent of parent of ....
  const Cell & ultimate_parent() const;
  // returns the parent of parent of parent of ....
  Cell & ultimate_parent();
  // return a vector of child cell indices
  inline const std::vector<size_t> & children() const { return m_children; }
  // returns a vector of indices of children of children of...
  std::vector<size_t> ultimate_children() const;

 protected:
  // this cell stuff
  std::size_t m_index;
  int m_marker, m_vtk_id;
  std::vector<std::size_t> m_vertices;
  std::vector<std::size_t> m_faces;
  // global grid stuff
  std::vector<Point> * pm_grid_vertices;
  std::vector<Cell> * pm_grid_cells;
  std::vector<Face> * pm_grid_faces;
  //  refinement
  std::size_t m_parent;
  std::vector<std::size_t> m_children;
  friend class Mesh;
  friend class active_cell_iterator;
  friend class active_cell_const_iterator;
};

}  // end namespace
