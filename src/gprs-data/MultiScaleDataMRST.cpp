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

void MultiScaleDataMRST::is2D_(std::vector<Point_3>& points)
{
	std::vector<double> x_,y_,z_;
	double rangex, rangey, rangez;

	for(auto pt:points)
	{
		x_.push_back(pt[0]);
		y_.push_back(pt[1]);
		z_.push_back(pt[2]);
	}

	std::vector<double>::iterator min, max;
	std::tie(min,max) = std::minmax_element(x_.begin(),x_.end());
	rangex = *max -*min;
	std::tie(min,max) = std::minmax_element(y_.begin(),y_.end());
	rangey = *max -*min;
	std::tie(min,max) = std::minmax_element(z_.begin(),z_.end());
	rangez = *max -*min;

	std::vector<double> ranges = {rangex,rangey,rangez};
	std::vector<std::size_t> sorted_ix = {0,1,2};
	//increasing inderected sort
	std::sort(sorted_ix.begin(),sorted_ix.end(),
			[&ranges](std::size_t i0, std::size_t i1){ return ranges[i0]>ranges[i1];});

	Vector_3 dl(0,0,0);
	switch(sorted_ix[2])
	{
	case 0:
		dl += Vector_3(ranges[sorted_ix[0]],0,0);
		break;
	case 1:
		dl += Vector_3(0,ranges[sorted_ix[0]],0);
		break;
	case 2:
		dl += Vector_3(0,0,ranges[sorted_ix[0]]);
		break;
	default:
		break;
	}

	std::cout << " Mesh can be considered 2D ... \n ranges: ";
		std::copy(ranges.begin(),ranges.end(),std::ostream_iterator<double>(std::cout," "));

		//TODO find a better criteria for 2D caracteriation
		//XXX work only if only one cells in z direction
		bool is2D = ranges[sorted_ix[2]] == 0 ;

	if(is2D)
	{

		std::cout << std::endl << ranges[sorted_ix[0]] << std::endl;
		std::cout << sorted_ix[2] << std::endl << dl << std::endl;

		int sz = x_.size();
		x_.insert(x_.end(),x_.begin(),x_.end());
		y_.insert(y_.end(),y_.begin(),y_.end());
		z_.insert(z_.end(),z_.begin(),z_.end());
		int i=0;

		switch(sorted_ix[2])
		{
		case 0:
			std::for_each(x_.begin(),x_.end(),[&sz,&i,&ranges,&sorted_ix](double& d){ d = (i++<sz) ? d+ranges[sorted_ix[0]] : d-ranges[sorted_ix[0]]; });
			break;
		case 1:
			std::for_each(y_.begin(),y_.end(),[&sz,&i,&ranges,&sorted_ix](double& d){ d = (i++<sz) ? d+ranges[sorted_ix[0]] : d-ranges[sorted_ix[0]]; });
			break;
		case 2:
			std::for_each(z_.begin(),z_.end(),[&sz,&i,&ranges,&sorted_ix](double& d){ d = (i++<sz) ? d+ranges[sorted_ix[0]] : d-ranges[sorted_ix[0]]; });
			break;

		default:
			break;
		}


		//reproduce points
		points.clear();
		for(int i=0; i<x_.size(); ++i)
			points.push_back( Point_3(x_[i],y_[i],z_[i]) );

	}

}

void MultiScaleDataMRST::convexHull(std::vector<Point_3>& points,const std::size_t block)
{
	CGAL::Object obj;
	is2D_(points);


	// compute convex hull of non-collinear points
	CGAL::convex_hull_3(points.begin(), points.end(), obj);
	if(const Polyhedron_3* poly = CGAL::object_cast<Polyhedron_3>(&obj) )
	{
		vHull_[block] = obj;
		//std::cout << "The convex hull contains " << poly->size_of_vertices() << " vertices" << std::endl;

		   //part for dbg
		      std::ofstream fout("part"+std::to_string(block)+".off",std::ofstream::out);
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
			std::cerr << "Planar hull" << std::endl;

		throw std::invalid_argument("Bad Value for Hull");

	}
}

void MultiScaleDataMRST::gatheredPoints()
{
	gathered_.resize(active_layer().cells_in_block.size());

	auto cMap = active_layer().block_internal_connections;
	for (auto it = cMap.begin(); it != cMap.end(); ++it)
	{

		const auto element = it.elements();

		//add all selected partition points
		for(auto elt_label:active_layer().cells_in_block[element.first])
			gathered_[element.first].push_back( toCGAL(grid.get_center(elt_label)));

		for(auto elt_label:active_layer().cells_in_block[element.second])
			gathered_[element.second].push_back( toCGAL(grid.get_center(elt_label)));

		// alt. use layer.coarse_to_fine
//		gathered_[element.first].push_back( toCGAL(grid.get_center( active_layer().coarse_to_fine[element.second] )) );
//		gathered_[element.second].push_back( toCGAL(grid.get_center( active_layer().coarse_to_fine[element.first] )) );
//		gathered_[element.first].push_back( toCGAL(grid.get_center( active_layer().coarse_to_fine[element.first] )) );
//		gathered_[element.second].push_back( toCGAL(grid.get_center( active_layer().coarse_to_fine[element.second] )) );
		// alt. using layer.block_centroids
		gathered_[element.first].push_back( toCGAL(active_layer().block_centroids[element.second]) );
		gathered_[element.second].push_back( toCGAL(active_layer().block_centroids[element.first]) );
		gathered_[element.first].push_back( toCGAL(active_layer().block_centroids[element.first]) );
		gathered_[element.second].push_back( toCGAL(active_layer().block_centroids[element.second]) );

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

// [deprecated] To check with my implemetation of AD-GPRS
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
