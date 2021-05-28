#pragma once

#include "ConnectionMap.hpp"
#include "algorithms/EdgeWeightedGraph.hpp"

#ifdef WITH_METIS
#include "metis.h"
#endif // WITH_METIS

#include <vector>
#include <exception>
#include <unordered_set>

namespace multiscale
{

using std::vector;
using std::size_t;

template <typename T>
using ConnectionMap = hash_algorithms::ConnectionMap<T>;


/* Wrapper class for METIS library */
class MetisInterface
{
 public:
  /* connections contains cell connections (or a graph)
   * n_blocks is how many partitions (coarse blocks) you wanna get
   * n_cells is how many unique untries are in the connection map;
   * if set to 0, then it's computed automatically */
  template <typename ConnType>
  static vector<size_t> build_partitioning(const ConnectionMap<ConnType> & connections,
                                           const size_t n_blocks,
                                           size_t n_cells = 0)
  {
#ifdef WITH_METIS
    if (n_blocks < 2) throw std::invalid_argument("number of blocks too small");
    if (n_cells == 0) n_cells = count_elements(connections);

    // generate input: xadj (similar to row_ptr in CSR format)
    vector<idx_t> xadj(n_cells + 1, 0);
    for (auto it = connections.begin(); it != connections.end(); ++it)
    {
        const auto elements = it.elements();
        const int ia = elements.first;
        const int ib = elements.second;
        ++xadj[ia+1];
        ++xadj[ib+1];
    }

    // convert to the same ascending format as row_ptr
    for(std::size_t ib = 0; ib < n_cells; ++ib)
      xadj[ib+1] += xadj[ib];

    // generate input: adj (similar to col_ind in CSR format)
    vector<idx_t> adj(xadj[n_cells]);
    // graph connection weights
    vector<idx_t> adj_weight(xadj[n_cells]);
    for (auto it = connections.begin(); it != connections.end(); ++it)
    {
        const auto elements = it.elements();
        const int ia = elements.first;
        const int ib = elements.second;
        adj_weight[xadj[ia]] = small_wgt; //small weight between blocks
        adj[xadj[ia]++] = ib;
        adj_weight[xadj[ib]] = small_wgt; //small weight between blocks
        // temporarily change the value of xadj for generating adj
        adj[xadj[ib]++] = ia;
    }

    //  restore xadj
    for(std::size_t ib = n_cells; ib != 0; --ib)
      xadj[ib] = xadj[ib-1];
    xadj[0] = 0;

    // partition the entire domain: see the user manual of METIS for details
    idx_t ncon = 1, objval = 0;
    idx_t options[METIS_NOPTIONS];
    METIS_SetDefaultOptions(options);
    options[METIS_OPTION_NCUTS] = 8;
    std::vector<idx_t> vwgt(n_cells, 1);
    std::vector<idx_t> size(n_cells, 1);
    idx_t icount = static_cast<idx_t>(n_cells);
    idx_t n_domains = static_cast<idx_t>(n_blocks);

    // output: the corresponding thread of each grid block (default: 0th thread)
    vector<idx_t> coarse_cell_id(n_cells, 0);
    // call METIS graph partitioner
    METIS_PartGraphKway(&icount, &ncon, &xadj[0], &adj[0], &vwgt[0], &size[0],
                        NULL/*&adjwgt[0]*/, &n_domains,
                        NULL, NULL, options, &objval,
                        &coarse_cell_id[0]);

    vector<size_t> partitioning(n_cells);
    for (std::size_t i=0; i<n_cells; ++i)
      partitioning[i] = static_cast<std::size_t>(coarse_cell_id[i]);

    return partitioning;

#else
  throw std::invalid_argument("METIS is not available, cannot perform partitioning");
#endif // WITH_METIS
  }  // end interface

  /*
  ** Partition graph g into n_blocks partitions.
  **
  ** The parameter imbalancing specifies the maximum allowed load imbalance among the partitions.
  ** A value of x indicates that the allowed load imbalance is (1 + x)/1000.
  ** The load imbalance for the jth constraint is defined to be
  ** maxi(w[j, i]) / t[j, i]), where w[j, i] is the fraction of the overall weight of the jth constraint that is as-
  ** signed to the ith partition and t[j, i] is the desired target weight of the jth constraint for the ith partition
  ** (i.e., that specified via -tpwgts)
  ** The value < 0 uses the default METIS imbalancing x = 30
   */
  static std::vector<size_t> partition(algorithms::EdgeWeightedGraph const & g, size_t n_blocks,
                                       double imbalancing = -1);

 private:
  MetisInterface() = delete;

  // computes number of unique elements in connection map
  template <typename ConnType>
  static size_t count_elements(const ConnectionMap<ConnType> & connections)
  {
    std::unordered_set<size_t> uniques;
    for (auto it = connections.begin(); it != connections.end(); ++it)
    {
      const auto elements = it.elements();
      uniques.insert(elements.first);
      uniques.insert(elements.second);
    }
    return uniques.size();
  }

  static const int small_wgt = 1;
  static const int large_wgt = 1000000;
};

}
