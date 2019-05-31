#pragma once

#include <iostream>
#include <fstream>
#include <vector>
#include <unordered_set>
#include <metis.h>
#include <assert.h>
#include "ConnectionMap.hpp"

using std::vector;
using std::cout;
using std::endl;




namespace mesh
{
// overload std::exception for METIS
class metis_error_input : public std::exception
{
	const char* what() const throw(){ return "Metis Bad Input" ; };
};
class metis_error_mem : public std::exception
{
	const char* what() const throw(){ return "Metis Bad Mem Usage" ;}
};
class metis_error : public std::exception
{
	const char* what() const throw(){ return "Metis Unknown Error" ;}
};

//METIS input data structure
class  MetisData
{

public:
	MetisData() = default;
	virtual ~MetisData(){};

	//exported from AD-GPRS and slightly modified
	void partitionConnectionList(const std::size_t           n_blocks/*num of partitions*/,
			const PureConnectionMap   & connection_list,
			std::vector<std::size_t>  & coarse_cell_idx,
			const std::size_t           n_elements = 0/*graph size*/);

	void export_METIS_partitions(const std::string& fname = "OUTPUT.METIS.txt");


private :

	// The number of balancing constraints. It should be at least 1
	idx_t ncon = 1, objval = 0;
	// The adjacency structure of the graph
	std::vector<idx_t> xadj, adj;
	// The weights of the vertices as described in Section 5.5.
	std::vector<idx_t> vwgt;
	// The size of the vertices for computing the total communication
	// volume as described in Section 5.7.
	std::vector<idx_t> size;

	/// What is actually in METIS
	idx_t* options_ = new idx_t[METIS_NOPTIONS];
	idx_t icount = 0;
	idx_t n_domains = 0;

	// consecutive weighting input/output args are left
	// unset by default for now
	//  vwgt (NULL) vsize (NULL) adjwgt (NULL)
	// tpwgts (NULL) ubvec (NULL)
	std::vector<std::size_t> coarse_cell_idx_;



	//helpers
	//set the number of vertexes
	MetisData& set_nvtxs(idx_t nvtxs) {
		icount = nvtxs;
		assert(icount); //reserve exception for more complex check
		return *this;
	};

	MetisData& set_ndom(idx_t ndom) {
		n_domains = ndom;
		assert(n_domains);
		return *this;
	};

	std::size_t count_elements(const PureConnectionMap& connection_list_) const;
	void process_CSRadjacency(const PureConnectionMap&);

}; // end of struct

}  // end namespace


