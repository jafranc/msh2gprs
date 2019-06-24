#pragma once

#include "Mesh.hpp"

//CGAL deps
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>

#include <vector>
#include <fstream>
#include <iostream>

//abstract class for support region
class SupportRegion{

public:
  SupportRegion(const Mesh* pm,
                const std::vector< std::vector<std::size_t> >& p):
    pMesh_(pm),partition_(p){};

protected:
  std::unique_ptr<Mesh> pMesh_;
  std::vector< std::vector<std::size_t> > partition_;

  typedef CGAL::Exact_predicates_inexact_constructions_kernel  K;
  typedef CGAL::Polyhedron_3<K>                     Polyhedron_3;
  typedef K::Point_3                                Point_3;
  typedef K::Segment_3                              Segment_3;
  typedef K::Triangle_3                             Triangle_3;

};

class MSRSBSupport : public SupportRegion
{
public :
  MSRSBSupport(const Mesh* pm,
               const std::vector< std::vector<std::size_t> >& p):
    SupportRegion(pm,p){};


private:
  std::vector< std::vector<Point_3> > gathered_;


  //converter (maybe unused)
  std::vector<Point_3> toCGAL(const angem::PointSet<3,double>& points){
    std::vector<Point_3> ptlist;
    for(auto it=points.begin(); it!=points.end(); ++it)
      ptlist.push_back(*it);

    return ptlist;
  };

  void getCentroid(const std::vector<std::size_t>& ix_list);

  void convexHull(const std::vector<Point_3>& points){

    CGAL::Object obj;
    // compute convex hull of non-collinear points
    CGAL::convex_hull_3(points.begin(), points.end(), obj);
    if(const Polyhedron_3* poly = CGAL::object_cast<Polyhedron_3>(&obj))
      std::cout << "The convex hull contains " << poly->size_of_vertices() << " vertices" << std::endl;

  };

  void gatheredPoints();
  // specific treatment at boundary
  void extrudeBoundaryFaceCenters();
  void edgeBoundaryCenters();

  //TODO specific geometric mean correction
  // void geometricMean();
};
