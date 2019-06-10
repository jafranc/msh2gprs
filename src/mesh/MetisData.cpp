/*
 * MetisData.cpp
 *
 *  Created on: May 29, 2019
 *      Author: jacques
 */

#include "MetisData.hpp"

namespace mesh {

//std::size_t MetisData::count_elements(const PureConnectionMap& connection_list) const
//{
//	std::unordered_set<std::size_t> set_elements;
//	for (auto it = connection_list.begin(); it != connection_list.end(); ++it)
//	{
//		const auto elements = it.elements();
//		set_elements.insert(elements.first);
//		set_elements.insert(elements.second);
//	}
//	return set_elements.size();
//
//	//TODO: replace by newly ConnectionData member function
//}

// generate input: xadj (similar to row_ptr in CSR format)
// generate input: adj (similar to col_ind in CSR format)
// TODO : try to improve const-ness
void MetisData::process_CSRadjacency(const PureConnectionMap& connection_list)
{
	xadj.resize(icount + 1);

	for (auto it = connection_list.begin(); it != connection_list.end(); ++it)
	{
		const auto elements = it.elements();
		const int ia = elements.first;
		const int ib = elements.second;
		xadj[ia+1]++;
		xadj[ib+1]++;
	}

	for(std::size_t ib = 0; ib < icount; ++ib)  // convert to the same ascending format as row_ptr
		xadj[ib+1] += xadj[ib];

	adj.resize(xadj[icount]);

	for (auto it = connection_list.begin(); it != connection_list.end(); ++it)
	{
		const auto elements = it.elements();
		const int ia = elements.first;
		const int ib = elements.second;
		adj[xadj[ia]++] = ib;
		// temporarily change the value of xadj for generating adj
		adj[xadj[ib]++] = ia;
	}

	// restore xadj
	for(std::size_t ib = icount; ib != 0; --ib)
		xadj[ib] = xadj[ib-1];
	xadj[0] = 0;

}
void MetisData::partitionConnectionList(const std::size_t           n_blocks,
		const PureConnectionMap & connection_list,
		const std::size_t           n_elements)
{
	// size of partitioning graph (icount) and number of partition (n_domains)
	set_nvtxs( static_cast<idx_t>((n_elements != 0) ? n_elements : connection_list.count_elements() ) );
	set_ndom( static_cast<idx_t>(n_blocks) );

	// generate input: xadj (similar to row_ptr in CSR format)
	process_CSRadjacency(connection_list);

	// partition the entire domain: see the user manual of METIS for details
	vwgt.resize(icount,1);
	size.resize(icount,1);

	// output: the corresponding thread of each grid block (default: 0th thread)
	vector<idx_t> coarse_cell_id(icount, 0);
	// call METIS graph partitioner

	try{
		// main triggers
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

	// copy to result vector
	coarse_cell_idx_.assign(coarse_cell_id.begin(),coarse_cell_id.end());

}; //end function

//writing METIS data partition formatted METIS like as
// a list of elt labels each lines standing for a different
// partition
void MetisData::write_METIS_partitions(const std::string& fname)
{
	//for debug purpose
		 int c=0;
		 std::vector< std::vector<std::size_t> > part(n_domains);
		 for(const auto& v : coarse_cell_idx_) part[v].push_back(c++);


		 std::ofstream fout(fname, std::ofstream::out);
     fout << "METISPART" << endl;
     for(const auto& p : part)
		 {
			 for(const auto& v : p) fout << v << " ";
			 fout << "/" << endl;
		 }
     fout << "/" << endl;
		 fout.close();
		 //end for debug purpose
}


}; // end of namespace
