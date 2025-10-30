#include "D:\Projects\C++\ns_fda\headers\inout.hpp"
#include "D:\Projects\C++\ns_fda\headers\mesh_n_model.hpp"
int main() {
  model_data model;
  const std::string inputFile = "D:\\Projects\\C++\\ns_fda\\config.txt";
  int code = input(inputFile, model);
  std::cout << "Domain partition: " << model.domainPartition << '\n';
  std::cout << "X_axis length: " << model.xLen << '\n';
  std::cout << "Y_axis length: " << model.yLen << '\n';
  std::cout << "Z_axis length: " << model.zLen << '\n';
  std::cout << "Reynild's number: " << model.Reyn << '\n';

  std::vector<double> xMesh;
  meshInit(xMesh, model, model.xLen);
  std::vector<double> yMesh;
  meshInit(yMesh, model, model.yLen);
  std::vector<double> zMesh;
  meshInit(zMesh, model, model.zLen);

  const std::string outputMeshFile = "D:\\Projects\\C++\\ns_fda\\";
  code = meshOutput(outputMeshFile + "xMesh.txt", xMesh);
  code = meshOutput(outputMeshFile + "yMesh.txt", yMesh);
  code = meshOutput(outputMeshFile + "zMesh.txt", zMesh);
}
//=================================================================================================