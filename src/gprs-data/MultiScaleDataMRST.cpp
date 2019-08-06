#include "MultiScaleDataMRST.hpp"

namespace multiscale {

using Point = angem::Point<3,double>;
using std::unordered_set;

MultiScaleDataMRST::MultiScaleDataMRST(mesh::Mesh  & grid, const size_t  n_blocks)
:
						   MultiScaleData(grid,n_blocks)
{
	vHull_.resize(n_blocks);
	gathered_.resize(n_blocks);
	active_layer().support_internal.resize(n_blocks);
	active_layer().support_overlap.resize(n_blocks);
	active_layer().support_boundary.resize(n_blocks);
}

// main method that identifies regions where shape functions exist
void MultiScaleDataMRST::build_support_regions()
{

	std::cout << "finding block centroids...";
	find_centroids();
	std::cout << "OK" << std::endl;

	std::cout << "building block connections...";
	build_block_connections();
	std::cout << "OK" << std::endl;


	std::cout << "build support regions..." << std::flush;

	for (std::size_t block = 0; block < active_layer().n_blocks; ++block)
	{
		std::cout << "\rbuild support regions... block "
				<< (int)(100*block/active_layer().n_blocks)  << "% " <<std::flush;
		build_support_region(block);
	}

	std::cout << "OK" << std::endl;

}

void MultiScaleDataMRST::build_support_region(const std::size_t block)
{
	//TODO
	if( !vHull_[block].empty() )
	{
		//add all cells in the parition
		support_region_impl_(block);
	}
	else
	{
		if(gathered_[block].empty())
			gatheredPoints();
		else
			convexHull(gathered_[block],block);

		build_support_region(block);
	}
}

void MultiScaleDataMRST::convexHull(const std::vector<Point_3>& points,const std::size_t block)
{
	CGAL::Object obj;
	// compute convex hull of non-collinear points
	CGAL::convex_hull_3(points.begin(), points.end(), obj);
	if(const Polyhedron_3* poly = CGAL::object_cast<Polyhedron_3>(&obj) )
		//  || const Triangle_3* t = CGAL::object_cast<Triangle_3>(&obj) ) //TODO implement 2D
	{
		vHull_[block] = obj;
		//std::cout << "The convex hull contains " << poly->size_of_vertices() << " vertices" << std::endl;
	}
	else
	{
		if(const Segment_3* p = CGAL::object_cast<Segment_3>(&obj))
			std::cerr << "Segment Hull " <<std::endl;
		else if(const Point_3* p = CGAL::object_cast<Point_3>(&obj))
			std::cerr << "Point Hull "<<std::endl;
		else if(const Triangle_3* p = CGAL::object_cast<Triangle_3>(&obj))
			std::cerr << "2D triangle Hull (to debug pending)" << std::endl; //TODO implement 2D

		throw std::invalid_argument("Bad Value for Hull");

	}
}

void MultiScaleDataMRST::gatheredPoints()
{
	gathered_.resize(active_layer().cells_in_block.size());

	auto cMap = active_layer().block_internal_connections;//grid.get_coarseConnectionMap();
	for (auto it = cMap.begin(); it != cMap.end(); ++it)
	{

		const auto element = it.elements();

		//add all selected partition points
		for(auto elt_label:active_layer().cells_in_block[element.first])
			gathered_[element.first].push_back( toCGAL(grid.get_center(elt_label)));

		for(auto elt_label:active_layer().cells_in_block[element.second])
			gathered_[element.second].push_back( toCGAL(grid.get_center(elt_label)));

		//cross add neighbors center
		gathered_[element.first].push_back( toCGAL(grid.get_center(active_layer().cells_in_block[element.second][0])) );
		gathered_[element.second].push_back( toCGAL(grid.get_center(active_layer().cells_in_block[element.first][0])) );

		//now boundaries (uncomment if needed)
		//gathered_[element.first].push_back( extrudeBoundaryFaceCenters(element.first) );
		//gathered_[element.second].push_back( extrudeBoundaryFaceCenters(element.second) );
	}

}


void MultiScaleDataMRST::fill_output_model(MultiScaleOutputData & model,
		const int layer_index) const
{
	const auto & layer = layers[layer_index];

	model.cell_data = true;
	model.n_coarse = layer.n_blocks;

	// partitioning
	model.partitioning.resize(layer.partitioning.size());
	std::copy(layer.partitioning.begin(), layer.partitioning.end(), model.partitioning.begin());

	//  centroids
	model.centroids = layer.coarse_to_fine;

	// support boundary
	model.support_boundary = layer.support_boundary;

	//  support internal
	model.support_internal = layer.support_internal;

	// support overlap
	model.support_overlap = layer.support_overlap;
}

void MultiScaleDataMRST::support_region_impl_(const std::size_t block)
{
	//CGAL returned lambda for the ray casting test
	CGAL::Side_of_triangle_mesh<Polyhedron_3,K> is_inside(*CGAL::object_cast<Polyhedron_3>(&vHull_[block]));

	auto& cMap = active_layer().block_internal_connections;

	vector<std::size_t>  vSupport;

	//start adding to the interior all the METIS partition points
	vSupport.insert(vSupport.begin(),
			active_layer().cells_in_block[block].begin(),active_layer().cells_in_block[block].end());
	//then testing if inside convex hull
	vector<std::size_t> neighs = cMap.get_neighbors(block);
	for(auto blockn : neighs)
		for(auto fc : active_layer().cells_in_block[blockn])
		{
			bool on_bounded = ( is_inside(toCGAL(grid.get_center(fc)) ) == CGAL::ON_BOUNDED_SIDE ||
					is_inside(toCGAL(grid.get_center(fc)) ) == CGAL::ON_BOUNDARY  );

			if(on_bounded) vSupport.push_back(fc);

		}//end for for loops

	assert( vSupport.size() > active_layer().cells_in_block[block].size() );

	//post treat support with non neighbors test
	for(auto fcI = vSupport.begin() ; fcI != vSupport.end(); ++fcI)
	{
		std::vector<std::size_t> mneigh = grid.get_neighbors(*fcI);
		for(auto neigh : mneigh)
		{//for all neigh if not found in subset then it is border
			//TODO try union find on this
			if( std::find(vSupport.begin(),vSupport.end(),neigh)
			== vSupport.end())
			{
				active_layer().support_boundary[block].insert(neigh);
				assert( active_layer().partitioning[neigh] != block );
			}
		}

	}


	//split interior from overlap
	for(auto it = vSupport.begin(); it != vSupport.end(); ++it)
	{
		if( active_layer().partitioning[*it]==block )
			active_layer().support_internal[block].insert(*it);
		else
			active_layer().support_overlap[block].insert(*it);
	}

}


//output operator for debug and temp solution export to gprs
void MultiScaleDataMRST::printFiles() const
{
	//output 3 needed files METIS.OUTPUT.txt MCONN.OUTPUT.txt MSRSB.OUTPUT.txt
	// METIS.OUTPUT.txt
	std::ofstream ofs("METIS.OUTPUT.txt");
	ofs << "METISPART" << std::endl;
	for(int block=0; block< active_layer().n_blocks; ++block)
		{
		std::copy(active_layer().support_internal[block].begin(),
						  active_layer().support_internal[block].end(),
						  std::ostream_iterator<std::size_t>(ofs," "));

		ofs<< " / " << std::endl;
		}
	ofs<< " / " << std::endl;
	ofs.close();

	//MCONN.OUTPUT.txt
	ofs.open("MCONN.OUTPUT.txt");
	ofs << "METISCONN" << std::endl;
	auto& cMap = active_layer().block_internal_connections;//grid.get_coarseConnectionMap();
	//check for all points in the neighbors
	for(int block=0; block< active_layer().n_blocks; ++block)
	{
		vector<std::size_t> neighs = cMap.get_neighbors(block);
		ofs << neighs.size() << " ";

		std::copy(neighs.begin(),
				  neighs.end(),
				  std::ostream_iterator<std::size_t>(ofs," "));

		ofs<< " / " << std::endl;



	}
	ofs<< " / " << std::endl;
	ofs.close();

	//MSRSB.OUTPUT.txt
	// output: Nbound Ninternal internal _labels overlap_labels boundary_labels
	ofs.open("MSRSB.OUTPUT.txt");
	ofs << "MSRSBSupport" << std::endl;
	for(int block=0; block< active_layer().n_blocks; ++block)
	{
		ofs << active_layer().support_boundary[block].size() << " "
			<< active_layer().support_internal[block].size() << " ";
			//<< active_layer().support_overlap[block].size() << " ";

		std::copy(active_layer().support_internal[block].begin(),
				  active_layer().support_internal[block].end(),
				  std::ostream_iterator<std::size_t>(ofs," "));

		std::copy(active_layer().support_overlap[block].begin(),
						  active_layer().support_overlap[block].end(),
						  std::ostream_iterator<std::size_t>(ofs," "));

		std::copy(active_layer().support_boundary[block].begin(),
						  active_layer().support_boundary[block].end(),
						  std::ostream_iterator<std::size_t>(ofs," "));
		ofs<< " / " << std::endl;
	}
	ofs<< " / " << std::endl;
	ofs.close();
}


}//end namespace multiscale
