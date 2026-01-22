#pragma once
#include <iostream>
#include <vector>
#include <cmath>
using uint = unsigned int;
// Modeling data
struct model_data {
  uint fdaNumber       {}; // number of finite difference approximation scheme
  uint domainPartition {}; // domain partition
  uint timePartition   {}; // time partition
  double xLen          {}; // domain x_axis length
  double yLen          {}; // domain y_axis length
  double zLen          {}; // domain z_axis length
  double duration      {}; // flowing duration
  double Reyn          {}; // Reynolds number

  const std::string PATH = "D:\\Projects\\C++\\ns_fda\\";
  const std::string PATH_log = "\\log.txt";
  const std::string PATH_residual = "\\residual.txt";
  const std::string PATH_config = "\\config";
  void operator=(const model_data&& other) 
  {
    fdaNumber = other.fdaNumber; 
    domainPartition = other.domainPartition; 
    timePartition = other.timePartition;
    xLen = other.xLen; yLen = other.yLen; zLen = other.zLen; 
    duration = other.duration; Reyn = other.Reyn;
  }
  void show() const 
  {
    std::cout << "X x Y x Z: " << xLen << " x " << yLen << " x " << zLen << '\n';
    std::cout << "flowing duration: " << duration << '\n';
    std::cout << "domain partition fineness: " << domainPartition << '\n';
    std::cout << "time partition fineness: " << timePartition << '\n';
    std::cout << "Reynolds number: " << Reyn << '\n';
  }

};
//=================================================================================================