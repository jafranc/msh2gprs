#pragma once

#include "mesh/Mesh.hpp"
#include "MetisInterface.hpp"
#include "LayerDataMSRSB.hpp"
#include "MultiScaleOutputData.hpp"
#include "UnionFindWrapper.hpp"
#include "tuple_hash.hpp"
#include <algorithm>  // std::max_element
#include <vector>

namespace multiscale
{

using std::size_t;
using std::vector;

// Abstract class for Multiscale Partitioning
// holds Metis related Mechanism
class MultiScaleData
{

public:
	MultiScaleData( mesh::Mesh  & grid, const size_t  n_blocks);
	virtual ~MultiScaleData(){};


	void build_data();

	//Talk-to-METIS methods
	// call to metis to obtain partitioning
	void build_partitioning();
	// build connection list and find centroids
	// build inverted partitioning block -> cells
	void build_cells_in_block();
	// build connection map that stores faces between blocks and their centers
	void build_block_connections();
	// find geometric centers of coarse blocks
	void find_centroids();

	// main method that identifies regions where shape functions exist
	virtual void build_support_regions() = 0;
	virtual void fill_output_model(MultiScaleOutputData & model, const int layer_index = 0) const = 0;



protected:
	// get reference to the active layer
	inline
	LayerDataMSRSB & active_layer(){return layers[active_layer_index];}
	const LayerDataMSRSB & active_layer() const {return layers[active_layer_index];}

	const mesh::Mesh & grid;
	vector<LayerDataMSRSB> layers;
	size_t active_layer_index;


private:


	virtual void build_support_region(const std::size_t block) = 0;


};

}//end namespace multiscale
