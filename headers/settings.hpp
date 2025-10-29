#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
// Modeling data
struct model_data {
  int domainPartition {}; // domain partition
  double xLen         {}; // domain x_axis length
  double yLen         {}; // domain y_axis length
  double zLen         {}; // domain z_axis length
  double Reyn         {}; // Reynolds number
};