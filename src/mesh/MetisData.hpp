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
  class metis_empty_error : public std::exception
  {
    const char* what() const throw(){ return "Metis Not runned Error" ;}
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
                                 const std::size_t           n_elements = 0/*graph size*/);

    //writing METIS data partition formatted METIS like as
    // a list of elt labels each lines standing for a different
    // partition
    // NB: wrap in mesh::Mesh for direct access
    void write_METIS_partitions(const std::string& fname);

    //getter
    const std::vector<std::size_t>& getCoarseCellIdx() const { return coarse_cell_idx_;};
    const std::vector< std::vector<std::size_t> >& getPart() const  { return part_;};
    idx_t get_ndom() const { return n_domains;};

    ///use to move front centroid cell label in output
    //TODO check memory sanity
    void move_front(std::size_t coarse_cell, std::size_t index )
    {
      auto found = std::find(part_[coarse_cell].begin(),part_[coarse_cell].end(),index);
      if(found != part_[coarse_cell].end())
        {
          part_[coarse_cell].insert(part_[coarse_cell].begin(),*found);
          part_[coarse_cell].erase(++found);
        }
    }

    // XXX: erase in Release
    void dbg_print_adj()
    {
      cout << " printing xadj :\n";
      for(const auto& v : xadj) cout << v << " ";
      cout <<endl;

      cout << " printing adj :\n";
      for(const auto& v : adj) cout << v << " ";
      cout <<endl;
    }

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
    std::vector<std::size_t> coarse_cell_idx_;//raw METIS output
    std::vector< std::vector<std::size_t> > part_;//coarse to fine mapping (needed for ordering)

    //helpers

    //setter the number of vertexes
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

    void process_CSRadjacency(const PureConnectionMap&);

  }; // end of struct

}  // end namespace
