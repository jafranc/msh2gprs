#include "SupportRegion.hpp"


//this is used to fill gathered_ list per partition
// this should include:
//  * faces and extrudedBoundary faces centroids
//  * blocks centroids and displaced boundary block centroids
// in order to hullify them

void MSRSBSupport::gatheredPoints()
{
  gathered_.resize(partition_.size());

      auto cMap = pMesh_->get_coarseConnectionMap();
      for (auto it = cMap->begin(); it != cMap->end(); ++it)
      {

        const auto element = it.elements();

        //add its center and neigh center to its list (resp. neigh list)
        for(auto elt_label:partition_[element.first])
          gathered_[element.first].push_back( toCGAL(pMesh_->get_center(elt_label)));


        for(auto elt_label:partition_[element.second])
          gathered_[element.second].push_back( toCGAL(pMesh_->get_center(elt_label)));

        //cross add
        gathered_[element.first].push_back( toCGAL(pMesh_->get_center(partition_[element.second][0])) );
        gathered_[element.second].push_back( toCGAL(pMesh_->get_center(partition_[element.first][0])) );

        //now boundaries
        gathered_[element.first].push_back( extrudeBoundaryFaceCenters(element.first) );

        gathered_[element.second].push_back( extrudeBoundaryFaceCenters(element.second) );
      }

}

MSRSBSupport::Point_3 MSRSBSupport::extrudeBoundaryFaceCenters(std::size_t I)
{
  auto cB = pMesh_->get_cBoundary();
    Point_3 cc = toCGAL( pMesh_->get_center(partition_[I][0]) );
  if(!cB.empty())
  {
    Point_3 fc = getCentroid(cB[I]);
    //centroids always strore first

  return fc+2*(fc-cc);
  }
  else
    return cc;
};

// void MSRSBSupport::edgeBoundaryCenters()
// {
// };

//get coarse face centroids
MSRSBSupport::Point_3 MSRSBSupport::getCentroid(const std::vector<std::size_t>& ix_list)
{
  angem::Point centroids;

  for(auto it=ix_list.begin(); it!=ix_list.end(); ++it)
    {
      centroids += pMesh_->get_center(*it)/ix_list.size();
    }

  return Point_3(centroids[0],centroids[1],centroids[2]);
}

void MSRSBSupport::convexHull(const std::vector<Point_3>& points, int offset_/* dbg only*/)
{

  CGAL::Object obj;
  // compute convex hull of non-collinear points
  CGAL::convex_hull_3(points.begin(), points.end(), obj);
  if(const Polyhedron_3* poly = CGAL::object_cast<Polyhedron_3>(&obj) )
    //  || const Triangle_3* t = CGAL::object_cast<Triangle_3>(&obj) ) // to debug
    {
      vHull_.push_back(obj);
      std::cout << "The convex hull contains " << poly->size_of_vertices() << " vertices" << std::endl;

      //part for dbg
      std::ofstream fout("part"+std::to_string(offset_)+".off",std::ofstream::out);
      fout << *poly << std::endl;
      fout.close();
    }
  else
    {
      if(const Segment_3* p = CGAL::object_cast<Segment_3>(&obj))
        std::cerr << "Segment Hull " <<std::endl;
      else if(const Point_3* p = CGAL::object_cast<Point_3>(&obj))
        std::cerr << "Point Hull "<<std::endl;
      else if(const Triangle_3* p = CGAL::object_cast<Triangle_3>(&obj))
        std::cerr << "2D triangle Hull (to debug pending)" << std::endl;

      throw std::invalid_argument("Bad Value for Hull");

    }
};

void MSRSBSupport::genSupport(std::size_t I)
{

  if( !vHull_[I].empty() )
    {
      //needed lamda to test
      CGAL::Side_of_triangle_mesh<Polyhedron_3,K> is_inside(*CGAL::object_cast<Polyhedron_3>(&vHull_[I]));
      //add all cells in the parition
      std::get<0>(vSupport_[I].first) = 0;
      std::get<1>(vSupport_[I].first) = 0;//partition_[I].size();
      vSupport_[I].second.insert( vSupport_[I].second.end(), partition_[I].begin(),partition_[I].end());

      auto cMap = pMesh_->get_coarseConnectionMap();
      //check for all points in the neighbors
      vector<std::size_t> neighs = cMap->get_neighbors(I);
      for(auto In : neighs)
        for(auto fc : partition_[In])
          {

             bool on_bounded = ( is_inside(toCGAL(pMesh_->get_center(fc)) ) == CGAL::ON_BOUNDED_SIDE || is_inside(toCGAL(pMesh_->get_center(fc)) ) == CGAL::ON_BOUNDARY );

              if(on_bounded) vSupport_[I].second.push_back(fc);

          }//end for for loops

      assert( vSupport_[I].second.size() > partition_[I].size() );

      //post treat suuport with non neighbors test
      std::vector<std::size_t> bc;
      for(auto fcI = vSupport_[I].second.begin() ; fcI != vSupport_[I].second.end(); ++fcI)
      {
        for(auto neigh:pMesh_->get_neighbors(*fcI))
          {//for all neigh if not found in subset then it is border
            //TODO try union find on this
            if( std::find(vSupport_[I].second.begin(),vSupport_[I].second.end(),neigh) == vSupport_[I].second.end())
              {
                bc.push_back(*fcI);
                vSupport_[I].second.erase(fcI);
                --fcI;
                break;
              }
            }

      }

      //can be restrained to I's neighbors
      auto fit = vSupport_[0].second.begin();
      for(int In=0; In < vSupport_.size(); In++ )
        {
          for(auto v:pMesh_->get_cBoundary()[In])
            if( (fit = std::find(vSupport_[I].second.begin(),vSupport_[I].second.end(),v))
                != vSupport_[I].second.end() )
              {
                bc.push_back(*fit);
                vSupport_[I].second.erase(fit);

              }
        }

      //sort and unique border as a 3-outisde neigh cell is added three times
      {
        std::cout << " bc :" << bc.size() << std::endl;
        std::sort(bc.begin(),bc.end());
        std::vector<std::size_t>::iterator it = std::unique(bc.begin(),bc.end());
        bc.erase(it,bc.end());
        std::cout << " bc :" << bc.size() << std::endl;
      }

      //reorder boundary last
      //TODO: check if boundary cells and inter-coarse blocks cells cannot duplicates

      //move back bc
      /*int nb = move_back(vSupport_[I].second,bc);

      std::get<0>(vSupport_[I].first) += nb;

      //adding current to list to spare some lines
      for(int In=0; In < vSupport_.size(); In++ )
        {
          int nb = move_back(vSupport_[I].second,pMesh_->get_cBoundary()[In]);
          std::get<0>(vSupport_[I].first) += nb;
          if(I==In)  std::get<1>(vSupport_[I].first) -=nb;

        }*/

      //finally sort internal with I members first
      // and paste bc at the end
      //

      std::vector< std::size_t > metisVec = pMesh_->getCoarseCellIdx();

      std::sort(vSupport_[I].second.begin(),vSupport_[I].second.end(),
                [&I,&metisVec] (std::size_t v1, std::size_t v2)
                { return (metisVec[v1]-I)< (metisVec[v2]-I); } // as std::size_t is unsigned should work
                );

      std::get<1>(vSupport_[I].first) =
        std::count_if(vSupport_[I].second.begin(),vSupport_[I].second.end(),
                    [&metisVec,&I] (std::size_t v)
                    { return metisVec[v]==I; }
                    );

      std::get<0>(vSupport_[I].first) = bc.size();
      vSupport_[I].second.insert(vSupport_[I].second.end(),bc.begin(),bc.end());
    }//if hull non empty
};

//return number of swapped elmts
std::size_t MSRSBSupport::move_back(std::vector<std::size_t>& mainVec, const std::vector<std::size_t>& swapped_list)
{
  std::size_t nb=0;
  std::vector<std::size_t> found_val;

  for(auto val : swapped_list)
    {
      auto found = std::find(mainVec.begin(),mainVec.end(),val);
      if(found!=mainVec.end())
        {
          mainVec.erase(found);
          ++nb;
          found_val.push_back(val);
        }
    }

  mainVec.insert(mainVec.end(),found_val.begin(),found_val.end());

  return nb;
};
