#include "D:\Projects\C++\ns_fda\headers\inout.hpp"
#include "D:\Projects\C++\ns_fda\headers\mesh_n_model.hpp"
#include "D:\Projects\C++\ns_fda\headers\compute_flow.hpp"

//=================================================================================================
int main() {
  model_data model;
  const std::string inputFile = "D:\\Projects\\C++\\ns_fda\\config";
  int code = dataInput(inputFile, model);
  model.show();

  // mesh init (is it necessary?)
  std::vector<double> xMesh;
  meshInit(xMesh, model, model.xLen);
  std::vector<double> yMesh;
  meshInit(yMesh, model, model.yLen);
  std::vector<double> zMesh;
  meshInit(zMesh, model, model.zLen);

  // function init
  std::vector<double> velX;
  initialConditions(velX, 0, model);
  std::vector<double> velY;
  initialConditions(velY, 1, model);
  std::vector<double> velZ;
  initialConditions(velZ, 2, model);
  std::vector<double> pressure;
  initialConditions(pressure, 3, model);

  // flow computation
  compute(model, velX, velY, velZ, pressure);

  // final res output
  const std::string outputMeshFile = "D:\\Projects\\C++\\ns_fda\\";
  code = meshOutput(outputMeshFile + "xMesh.txt", xMesh);
  code = meshOutput(outputMeshFile + "yMesh.txt", yMesh);
  code = meshOutput(outputMeshFile + "zMesh.txt", zMesh);

  const std::string outputFuncFile = "D:\\Projects\\C++\\ns_fda\\";
  code = funcOutput(outputFuncFile, "v1", ".txt", velX, model);
  code = funcOutput(outputFuncFile, "v2", ".txt", velY, model);
  code = funcOutput(outputFuncFile, "v3", ".txt", velZ, model);
  code = funcOutput(outputFuncFile, "p", ".txt", pressure, model);

  std::ofstream outputFile("D:\\Projects\\C++\\ns_fda\\log.txt", std::ios::app);
  if (outputFile.is_open()) outputFile << "execution is finished" << '\n';
}
//=================================================================================================