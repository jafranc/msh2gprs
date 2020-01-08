#pragma once

#include "angem/Point.hpp"
#include <vector>

namespace discretization {

enum ConnectionType
{
  matrix_matrix = 1,
  matrix_fracture = 2,
  fracture_fracture = 3
};

struct ConnectionData
{
  ConnectionType type;
  std::vector<double> coefficients;
  std::vector<size_t> elements;
  double area;
  angem::Point<3,double> normal;
  angem::Point<3,double> center;
};

}  // end namespace discretization
