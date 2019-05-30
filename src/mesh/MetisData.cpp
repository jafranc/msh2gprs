/*
 * MetisData.cpp
 *
 *  Created on: May 29, 2019
 *      Author: jacques
 */

#include "MetisData.hpp"

namespace mesh {

void MetisData::partitionConnectionList(const std::size_t           n_blocks,
		const PureConnectionMap   & connection_list,
		std::vector<std::size_t>  & coarse_cell_idx,
		const std::size_t           n_elements)
{
	// size of partitioning graph
	std::size_t n_graph_nodes = n_elements;
	if (n_graph_nodes == 0)  // count nodes
	{
		std::unordered_set<std::size_t> set_elements;
		for (auto it = connection_list.begin(); it != connection_list.end(); ++it)
		{
			const auto elements = it.elements();
			set_elements.insert(elements.first);
			set_elements.insert(elements.second);
		}
		n_graph_nodes = set_elements.size();
	}

	// generate input: xadj (similar to row_ptr in CSR format)
	//TODO: transfer to METIS class
	xadj.resize(n_graph_nodes + 1);

	for (auto it = connection_list.begin(); it != connection_list.end(); ++it)
	{
		const auto elements = it.elements();
		const int ia = elements.first;
		const int ib = elements.second;
		xadj[ia+1]++;
		xadj[ib+1]++;
	}

	for(std::size_t ib = 0; ib < n_graph_nodes; ++ib)  // convert to the same ascending format as row_ptr
		xadj[ib+1] += xadj[ib];

	// generate input: adj (similar to col_ind in CSR format)
	//TODO: transfer to METIS class
	//TODO: at for bound check should be erased in release
	adj.resize(xadj.at(n_graph_nodes));

	for (auto it = connection_list.begin(); it != connection_list.end(); ++it)
	{
		const auto elements = it.elements();
		const int ia = elements.first;
		const int ib = elements.second;
		adj[xadj[ia]++] = ib;
		// temporarily change the value of xadj for genrating adj
		adj[xadj[ib]++] = ia;
	}

	// restore xadj
	for(std::size_t ib = n_graph_nodes; ib != 0; --ib)
		xadj[ib] = xadj[ib-1];
	xadj[0] = 0;

	// partition the entire  if(METIS_ERROR_INPUT) throw metis_error_input();
	//	  if(METIS_ERROR_MEMORY)throw metis_error_mem();
	//	  if(METIS_ERROR) throw metis_error(); domain: see the user manual of METIS for details
	vwgt.resize(n_graph_nodes,1);
	size.resize(n_graph_nodes,1);

	set_nvtxs( static_cast<idx_t>(n_graph_nodes) );
	set_ndom( static_cast<idx_t>(n_blocks) );

	// output: the corresponding thread of each grid block (default: 0th thread)
	vector<idx_t> coarse_cell_id(n_graph_nodes, 0);
	// call METIS graph partitioner

	try{
		// main triggerx
		METIS_SetDefaultOptions(options_);
		options_[METIS_OPTION_NCUTS] = 8;
		//// from METIS doc
		//	  METIS_API(int) METIS_PartGraphKway(idx_t *nvtxs, idx_t *ncon,
		//	  idx_t *xadj, idx_t *adjncy, idx_t *vwgt, idx_t *vsize,
		//	  idx_t *adjwgt, idx_t *nparts,
		//	  real_t *tpwgts, real_t *ubvec, idx_t *options, idx_t *edgecut,
		//	  idx_t *part);

		idx_t ret = METIS_PartGraphKway(&icount, &ncon,
				&xadj[0], &adj[0], &vwgt[0], &size[0],
				NULL, &n_domains,
				NULL, NULL, options_, &objval,
				&coarse_cell_id[0]);

		if(ret==METIS_ERROR_INPUT) throw metis_error_input();
		if(ret==METIS_ERROR_MEMORY)throw metis_error_mem();
		if(ret==METIS_ERROR) throw metis_error();

	}
	catch(const std::exception& e)
	{
		// METIS has actually error code too bad not to use it
		cout << "An error occur using METIS :" << e.what() << endl;
	}

	// copy to result vector //TODO overuse of two array, overload explicit conversion function
	coarse_cell_idx.resize(n_graph_nodes, 0);
	int i=0;
	std::for_each(coarse_cell_idx.begin(),coarse_cell_idx.end(),
			[&coarse_cell_id,&i](std::size_t& szi){ szi = static_cast<std::size_t>(coarse_cell_id[i++]); });
}; //end function



}; // end of namespace
