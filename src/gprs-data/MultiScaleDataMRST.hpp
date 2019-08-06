#pragma once

#include "mesh/Mesh.hpp"
#include "MultiScaleData.hpp"
#include "LayerDataMSRSB.hpp"
#include "MultiScaleOutputData.hpp"
#include "UnionFindWrapper.hpp"
#include "tuple_hash.hpp"
#include <algorithm>  // std::max_element
#include <vector>
#include <fstream>

//CGAL deps
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
//ref to CGAL::PMP/point_inside_example
#include <CGAL/Side_of_triangle_mesh.h>


namespace multiscale
{

using std::size_t;
using std::vector;

class MultiScaleDataMRST : public MultiScaleData
{
 public:
  MultiScaleDataMRST(mesh::Mesh  & grid, const size_t  n_blocks);
  virtual ~MultiScaleDataMRST(){};

  // main method that identifies regions where shape functions exist
  void build_support_regions();
  virtual void fill_output_model(MultiScaleOutputData & model, const int layer_index = 0) const;

  //gathered points of interest for convex hull computation
  virtual void gatheredPoints();

  void printFiles() const;

  private:

  typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
  typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
  typedef K::Point_3                                Point_3;
  typedef K::Segment_3                              Segment_3;
  typedef K::Triangle_3                             Triangle_3;

  std::vector< std::vector<Point_3> > gathered_;
  //CGAL::object can be either Triangle_3 or Polyhedron_3
    // for 2D and 3D case is else (Popint_3 or Segment_3)
    // throw a bad setup
    std::vector< CGAL::Object> vHull_;


  // build support region for a block
  void build_support_region(const std::size_t block);
  // Compute convexHull object from list of points thanks to CGAL implementation
  //and store it at vHull_[block]
  //TODO reimplenent 3D-Hull (apart from CGAL)
  void convexHull(const std::vector<Point_3>& points,const std::size_t block);
  //converter
  //TODO if CGAL kept turn into implicit conversion op
  inline  Point_3 toCGAL(const angem::Point<3,double>& pt)
  {
    return Point_3(pt[0],pt[1],pt[2]);
  };

  //seperated implementation for build support region
  void support_region_impl_(const std::size_t block);

};

}// end namespace multiscale
