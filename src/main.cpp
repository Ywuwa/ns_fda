#include "D:\Projects\C++\ns_fda\headers\inout.hpp"
int main() {
  model_data model;
  const std::string inputFile = "D:\\Projects\\C++\\ns_fda\\config.txt";
  int code = input(inputFile, model);
  std::cout << model.mesh_knots_num << '\n';
  std::cout << model.Reyn << '\n';
}
//=================================================================================================