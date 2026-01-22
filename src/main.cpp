#include "../headers/inout.hpp"
#include "../headers/mesh_n_model.hpp"
#include "../headers/compute_flow.hpp"
#include "../headers/concurrent_compute_flow.hpp"
#include <chrono>
#include <iomanip> // std::scientific, std::setprecision

//=================================================================================================
int main() {
  model_data model;
  const std::string inputFile = model.PATH + model.PATH_config;
  int code = dataInput(inputFile, model);
  model.show();

  // function init
  const uint vecSize ( (model.domainPartition + 1) * 
                        (model.domainPartition + 1) * 
                        (model.domainPartition + 1) );
  std::vector<double> velX(vecSize);
  initialConditions(velX, 0, model);
  std::vector<double> velY(vecSize);
  initialConditions(velY, 1, model);
  std::vector<double> velZ(vecSize);
  initialConditions(velZ, 2, model);
  std::vector<double> pressure(vecSize);
  initialConditions(pressure, 3, model);

  const std::string outputFuncFile = model.PATH;
  code = funcOutput(outputFuncFile, "v1", std::to_string(0), ".txt", velX, model, true);
  code = funcOutput(outputFuncFile, "v2", std::to_string(0), ".txt", velY, model, true);
  code = funcOutput(outputFuncFile, "v3", std::to_string(0), ".txt", velZ, model, true);
  code = funcOutput(outputFuncFile, "p", std::to_string(0), ".txt", pressure, model, true);

  auto start = std::chrono::high_resolution_clock::now();
  // flow computation
  switch (model.fdaNumber)
  {
  case 1:
    compute_cube_FDA1_3(model, velX, velY, velZ, pressure);
    break;
  case 2:
    compute_cube_FDA2_4(model, velX, velY, velZ, pressure);
    break;
  case 3:
    compute_cube_FDA1_3(model, velX, velY, velZ, pressure);
    break;
  case 4:
    compute_cube_FDA2_4(model, velX, velY, velZ, pressure);
    break;
  default:
    break;
  }
  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  std::cout << "Время выполнения: " << elapsed.count() << " секунд" << std::endl;

  // flow computation with OpenMP
  /*auto startConc = std::chrono::high_resolution_clock::now();
  //conc_compute_cube_FDA1_3(model, velX, velY, velZ, pressure);
  //conc_compute_cube_FDA2_4(model, velX, velY, velZ, pressure);

  //compute_precise(model, velX, velY, velZ, pressure);
  auto endConc = std::chrono::high_resolution_clock::now();

  std::chrono::duration<double> elapsedConc = endConc - startConc;
  std::cout << "Время выполнения: " << elapsedConc.count() << " секунд" << std::endl;
  std::cout << std::scientific << std::setprecision(6)
              << "aka " << elapsedConc.count() << " секунд\n";*/
  
  std::ofstream outputFile(model.PATH + model.PATH_log, std::ios::app);
  if (outputFile.is_open()) outputFile << "execution is finished" << '\n';
}
//=================================================================================================