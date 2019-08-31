#pragma once

#include "angem/Point.hpp"
#include "angem/Polygon.hpp"
#include <Face.hpp>
// #include <mesh_methods.hpp>
#include <vector>

namespace mesh
{

using Point = angem::Point<3,double>;
// using FaceMap = std::unordered_map<hash_type, Face>;

class face_iterator : public std::iterator<std::random_access_iterator_tag, Face>
{
 public:
  // Default constructor
  face_iterator(const std::size_t                     face_index,
                std::vector<Face>                   & grid_faces,
                std::vector<angem::Point<3,double>> & grid_vertices);
  // Copy constructor
  face_iterator(const face_iterator & other);

  // access operator
  Face * operator*() { return &(grid_faces[face_index]); }
  // comparison
  // returns true if both iterators point to the same cell
  bool operator==(const face_iterator & other) const;
  // returns true if the iterators point to different cells
  bool operator!=(const face_iterator & other) const;

  // SETTERS
  // assignment operator
  face_iterator & operator=(const face_iterator & other);

  // GETTERS
  // get center of mass of a face
  Point center() const;
  // get face normal vector
  Point normal() const;
  // get face marker, (-1) if not defined
  int marker();
  // get face index
  std::size_t index() const ;//{return face_it->second.index;}
  // get an index of the parent (master) face that existed before the split
  // same as index() if the face has not been split
  std::size_t master_index() const {return operator*()->master_face_index;}
  // get vtk id of the face
  int vtk_id() const ;//{return face_it->second.vtk_id;}
  // get vector of neighbor cell indices
  inline const std::vector<std::size_t> & neighbors() const {return operator*()->neighbor_cells;}
  // get vector of face vertex coordinates
  std::vector<Point> vertex_coordinates() const;
  // get vector of face vertex indices
  std::vector<std::size_t> vertex_indices() const;
  // get angem::Polygon from vertices
  angem::Polygon<double> polygon() const;
  // incrementing
  // increment operator
  face_iterator & operator++();
  // decrement operator
  face_iterator & operator--();

 private:
  const Face * operator*() const { return &(grid_faces[face_index]); }

 private:
  std::size_t                             face_index;
  std::vector<Face>                     & grid_faces;
  std::vector<angem::Point<3,double>>   & grid_vertices;
};

}
