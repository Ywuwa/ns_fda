#include "D:\Projects\C++\ns_fda\headers\inout.hpp"
#include "D:\Projects\C++\ns_fda\headers\mesh_n_model.hpp"
//=================================================================================================
int main() {
  model_data model;
  const std::string inputFile = "D:\\Projects\\C++\\ns_fda\\config";
  int code = dataInput(inputFile, model);
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

  std::vector<double> velX;
  initialConditions(velX, 0, model);
  std::vector<double> velY;
  initialConditions(velY, 1, model);
  std::vector<double> velZ;
  initialConditions(velZ, 2, model);
  std::vector<double> pressure;
  initialConditions(pressure, 3, model);

  const std::string outputMeshFile = "D:\\Projects\\C++\\ns_fda\\";
  code = meshOutput(outputMeshFile + "xMesh.txt", xMesh);
  code = meshOutput(outputMeshFile + "yMesh.txt", yMesh);
  code = meshOutput(outputMeshFile + "zMesh.txt", zMesh);

  const std::string outputFuncFile = "D:\\Projects\\C++\\ns_fda\\";
  code = funcOutput(outputFuncFile, "v1", ".txt", velX, model);
  code = funcOutput(outputFuncFile, "v2", ".txt", velY, model);
  code = funcOutput(outputFuncFile, "v3", ".txt", velZ, model);
  code = funcOutput(outputFuncFile, "p", ".txt", pressure, model);
}
//=================================================================================================