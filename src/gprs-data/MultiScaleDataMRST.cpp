#include "MultiScaleDataMRST.hpp"

namespace multiscale {

using Point = angem::Point<3,double>;
using std::unordered_set;

MultiScaleDataMRST::MultiScaleDataMRST(mesh::Mesh  & grid, const size_t  n_blocks)
 :
   MultiScaleData(grid,n_blocks)
   {}

 // main method that identifies regions where shape functions exist
 void MultiScaleDataMRST::build_support_regions()
 {

	  std::cout << "finding block centroids...";
	  find_centroids();
	  std::cout << "OK" << std::endl;

	  std::cout << "building block connections...";
	  build_block_connections();
	  std::cout << "OK" << std::endl;

	  /*build_cells_in_block();
	 build_block_connections();
	 find_centroids();*/

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


 void build_support_region(const std::size_t block)
 {
	 //TODO
 }

}//end namespace multiscale
