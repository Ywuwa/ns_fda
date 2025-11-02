#pragma once
#include <vector>
using uint = unsigned int;
// Modeling data
struct model_data {
  uint domainPartition {}; // domain partition
  double xLen                  {}; // domain x_axis length
  double yLen                  {}; // domain y_axis length
  double zLen                  {}; // domain z_axis length
  double Reyn                  {}; // Reynolds number
};