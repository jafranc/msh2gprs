#include "MultiScaleData.hpp"

namespace multiscale {

using Point = angem::Point<3,double>;
using std::unordered_set;


MultiScaleData::MultiScaleData(mesh::Mesh  & grid, const size_t  n_blocks)
:
				grid(grid),
				active_layer_index(0)
{
	auto & layer = layers.emplace_back();
	layer.index = 0;
	layer.n_blocks = n_blocks;
	layer.n_cells = grid.n_cells();
}

void MultiScaleData::build_data()
{
  std::cout << "building METIS partitioning...";
  build_partitioning();
  std::cout << "OK" << std::endl;

  build_cells_in_block();

  build_support_regions();
}

// call to metis to obtain partitioning
void MultiScaleData::build_partitioning()
{

	auto & layer = active_layer();
	PureConnectionMap cell_connections;

	for (auto it = grid.begin_faces(); it != grid.end_faces(); ++it)
	{
		const auto & neighbors = it.neighbors();
		if (neighbors.size() == 2)  // not a boundary face
			cell_connections.insert_connection( neighbors[0], neighbors[1] );
	}

	layer.partitioning = multiscale::MetisInterface<hash_algorithms::empty>
	::build_partitioning(cell_connections, layer.n_blocks, layer.n_cells);



}

// build inverted partitioning block -> cells
void MultiScaleData::build_cells_in_block()
{

	auto & layer = active_layer();
	layer.cells_in_block.resize(layer.n_blocks);
	for (std::size_t cell=0; cell<layer.n_cells; ++cell)
		layer.cells_in_block[layer.partitioning[cell]].push_back(cell);

}

void MultiScaleData::build_block_connections()
{
	auto & layer = active_layer();

	for (auto face = grid.begin_faces(); face != grid.end_faces(); ++face)
	{
		const auto & neighbors = face.neighbors();
		if (neighbors.size() == 2)  // not a boundary face
		{
			const std::size_t i1 = layer.partitioning[neighbors[0]];
			const std::size_t i2 = layer.partitioning[neighbors[1]];;
			if (i1 != i2)
				if (!layer.block_internal_connections.connection_exists(i1, i2))
					layer.block_internal_connections.insert_connection(i1, i2);
		}
	}
}

void MultiScaleData::find_centroids()
{
	auto & layer = active_layer();
	layer.block_centroids.resize(layer.n_blocks);
	layer.coarse_to_fine.resize(layer.n_blocks);

	if(!layer.cells_in_block.empty())
	{
		for(auto block_it = layer.cells_in_block.begin(); block_it != layer.cells_in_block.end(); ++block_it)
		{
			Point centroids;
			std::size_t block_ix = std::distance(layer.cells_in_block.begin(),block_it);
			for(auto cell = block_it->begin(); cell != block_it->end(); ++cell )
				layer.block_centroids[block_ix] += grid.get_center(*cell)/block_it->size();

			//store centroids in layer
			layer.coarse_to_fine[block_ix] = grid.findNearest(layer.block_centroids[block_ix]);

		}

	}


}


}//end namespace multiscale
