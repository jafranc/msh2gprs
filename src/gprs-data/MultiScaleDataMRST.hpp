#pragma once

#include "mesh/Mesh.hpp"
#include "MultiScaleData.hpp"
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

class MultiScaleDataMRST : public MultiScaleData
{
 public:
  MultiScaleDataMRST(mesh::Mesh  & grid, const size_t  n_blocks);
  virtual ~MultiScaleDataMRST(){};

  // main method that identifies regions where shape functions exist
  void build_support_regions();
  virtual void fill_output_model(MultiScaleOutputData & model, const int layer_index = 0) const;

  private:

  // build support region for a block
  void build_support_region(const std::size_t block);




};

  }// end namespace multiscale
