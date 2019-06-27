#pragma once

#include "Mesh.hpp"

//CGAL deps
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/convex_hull_3.h>
//ref to CGAL::PMP/point_inside_example
#include <CGAL/Side_of_triangle_mesh.h>

#include <vector>
#include <fstream>
#include <iostream>

//abstract class for support region
class SupportRegion{

public:
  SupportRegion(const mesh::Mesh& mesh_m,
                const std::vector< std::vector<std::size_t> >& p):
    partition_(p)
  {
    pMesh_ = std::make_shared<mesh::Mesh>(mesh_m);
  };

  virtual ~SupportRegion(){};

protected:
  std::shared_ptr<mesh::Mesh> pMesh_;
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
  // TODO uplift if working to Parent Class
  typedef std::vector< std::pair<std::size_t , std::vector<std::size_t> > > Support_container;


  MSRSBSupport(const mesh::Mesh& mesh_m,
               const std::vector< std::vector<std::size_t> >& p):
    SupportRegion(mesh_m,p){

    gatheredPoints();
    vSupport_.resize(p.size());
  };

  virtual ~MSRSBSupport(){};

  std::vector< std::vector<Point_3> > getGathered(){ return gathered_;};
  void process_supports()
  {
    for(int i=0; i<gathered_.size(); i++)
      {
        convexHull(gathered_[i], i);
        genSupport(i);
      }
  };

   friend std::ostream& operator<<(std::ostream& os, const MSRSBSupport& support)
   {
     for(const auto& sp: support.vSupport_ )
       {
        os << sp.first << " ";
        for(auto v:sp.second)
          os << v << " ";
        os<<std::endl;
      }

     return os;
   }

private:
  std::vector< std::vector<Point_3> > gathered_;
  //CGAL::object can be either Triangle_3 or Polyhedron_3
  // for 2D and 3D case is else (Popint_3 or Segment_3)
  // throw a bad setup
  std::vector< CGAL::Object> vHull_;

  // TODO uplift if working to Parent Class
  Support_container vSupport_;


  //converter (maybe unused)
  std::vector<Point_3> toCGAL(const angem::PointSet<3,double>& points){
    std::vector<Point_3> ptlist;
    for(auto i=0; i<points.size(); i++)
      ptlist.push_back({points[i][0],points[i][1],points[i][2]});

    return ptlist;
  };

  //converter (maybe unused)
  std::vector<Point_3> toCGAL(const std::vector< angem::Point<3,double> >& points){
    std::vector<Point_3> ptlist;
    for(auto i=0; i<points.size(); i++)
      ptlist.push_back(toCGAL(points[i]));

    return ptlist;
  };

  Point_3 toCGAL(const angem::Point<3,double>& pt){
    return Point_3(pt[0],pt[1],pt[2]);}

  Point_3 getCentroid(const std::vector<std::size_t>& ix_list);

  void convexHull(const std::vector<Point_3>& points, int offset_);
  void genSupport(std::size_t I);
  void gatheredPoints();

  // specific treatment at boundary
  // (cf. MRST implementation)
  Point_3 extrudeBoundaryFaceCenters(std::size_t I);
  //void edgeBoundaryCenters();

  //TODO specific geometric mean correction
  // void geometricMean();

  //utilities
  void move_back(std::vector<std::size_t>& mainVec, const std::vector<std::size_t>& swapped_list);


};
