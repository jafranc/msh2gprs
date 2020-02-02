#pragma once
#include "angem/Point.hpp"
#include <vector>

namespace gprs_data {

using Point = angem::Point<3,double>;

class FeValues
{
 public:
  // Constructor
  FeValues(const int element_type, const size_t n_elements);
  /* update function gradients and jacobians */
  void update(const size_t element_tag);
  // get value of shape i in integration point q
  double value(const size_t i, const size_t q) const;
  // get gradient of shape i in integration point q
  Point grad(const size_t i, const size_t q) const;
  // get determinant
  double JxW(const size_t q) const;
  // number of integration points
  size_t n_q_points() const;
  // number of vertices in the reference element
  size_t n_vertices() const;

 protected:
  void initialize_();
  void compute_shape_grads_();

 private:
  const int _element_type;
  const size_t _n_elements;
  int _n_comp;
  std::vector<double> _ref_points;   // integration on reference element
  std::vector<double> _weights;      // integration weights
  std::vector<double> _ref_gradients;     // basis function gradients [dxi_duj] on ref element
  // jacobians
  // [e1g1Jxu, e1g1Jyu, e1g1Jzu, e1g1Jxv, ..., e1g1Jzw, e1g2Jxu, ..., e1gGJzw, e2g1Jxu, ...]
  std::vector<double> _jacobians;
  std::vector<double> _determinants;
  std::vector<double> _true_points;  // integration points on real element
  // for only a single element at a time
  std::vector<double> _grad;  // shape gradients on real element
  std::vector<double> _inv_determinants;
};

}  // end namespace gprs_data
