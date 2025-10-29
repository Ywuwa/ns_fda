#include "D:\Projects\C++\ns_fda\headers\inout.hpp"
int main() {
  model_data model;
  const std::string inputFile = "D:\\Projects\\C++\\ns_fda\\config.txt";
  int code = input(inputFile, model);
  std::cout << "Domain partition: " << model.domainPartition << '\n';
  std::cout << "X_axis length: " << model.xLen << '\n';
  std::cout << "Y_axis length: " << model.yLen << '\n';
  std::cout << "Z_axis length: " << model.zLen << '\n';
  std::cout << "Reynild's number: " << model.Reyn << '\n';
}
//=================================================================================================