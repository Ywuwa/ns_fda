#include "../headers/inout.hpp"
#include "../headers/mesh_n_model.hpp"
#include "../headers/compute_flow.hpp"

//=================================================================================================
int main() {
  model_data model;
  const std::string inputFile = model.PATH + "\\config";
  int code = dataInput(inputFile, model);
  model.show();

  // function init
  std::vector<double> velX;
  initialConditions(velX, 0, model);
  std::vector<double> velY;
  initialConditions(velY, 1, model);
  std::vector<double> velZ;
  initialConditions(velZ, 2, model);
  std::vector<double> pressure;
  initialConditions(pressure, 3, model);

  const std::string outputFuncFile = model.PATH;
  code = funcOutput(outputFuncFile, "v1", std::to_string(0), ".txt", velX, model, true);
  code = funcOutput(outputFuncFile, "v2", std::to_string(0), ".txt", velY, model, true);
  code = funcOutput(outputFuncFile, "v3", std::to_string(0), ".txt", velZ, model, true);
  code = funcOutput(outputFuncFile, "p", std::to_string(0), ".txt", pressure, model, true);
  // flow computation
  compute(model, velX, velY, velZ, pressure);
  //compute_precise(model, velX, velY, velZ, pressure);
  
  std::ofstream outputFile(model.PATH + "\\log.txt", std::ios::app);
  if (outputFile.is_open()) outputFile << "execution is finished" << '\n';
}
//=================================================================================================