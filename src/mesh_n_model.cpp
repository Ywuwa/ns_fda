#include "D:\Projects\C++\ns_fda\headers\mesh_n_model.hpp"

//=================================== INITIAL FUCNTIONS ===========================================
using abc = ABC_Flow;
double abc::initV1(const double x, const double y, const double z)
{
  // u = A*sin(kz) + C*cos(ky)  
  return abc::A * std::sin(abc::k * z) + abc::C * std::cos(abc::k * y);
}
double abc::initV2(const double x, const double y, const double z)
{
  // v = B*sin(kx) + A*cos(kz)
  return abc::B * std::sin(abc::k * x) + abc::A * std::cos(abc::k * z);
}
double abc::initV3(const double x, const double y, const double z)
{
  // w = C*sin(ky) + B*cos(kx)
  return abc::C * std::sin(abc::k * y) + abc::B * std::cos(abc::k * x);
}
double abc::initPress(const double x, const double y, const double z)
{
  const double u = abc::initV1(x, y, z);
  const double v = abc::initV2(x, y, z);
  const double w = abc::initV3(x, y, z);
  // p = p0 - V^2 / 2 - P
  return abc::p0 - 0.5 * (u*u + v*v + w*w) - abc::P;
}
//=================================================================================================

//===================================== AXIS MESH INIT ============================================
void meshInit(std::vector<double>& axisMesh, const model_data& params, const double axisLen)
{
  const uint domainPartition = params.domainPartition;
  const double axisStep (axisLen / domainPartition);
  axisMesh.reserve(domainPartition + 1);
  for (auto i = 0; i < domainPartition + 1; i++)
  {
    axisMesh.emplace_back(i * axisStep);
  }
}
//=================================================================================================

//=================================== INITIAL CONDITIONS ==========================================
void initialConditions(
  std::vector<double>& feature, const uint initFuncInd, const model_data& params)
{
  const uint domainPartition = params.domainPartition;
  const double xStep (params.xLen / domainPartition);
  const double yStep (params.yLen / domainPartition);
  const double zStep (params.zLen / domainPartition);
  feature.reserve( (domainPartition + 1) * (domainPartition + 1) * (domainPartition + 1) );
  ABC_Flow functionSet;
  functionContainer fC(functionSet);

  for (auto k = 0; k < domainPartition + 1; k++)      // Z-Axis
  {
    for (auto j = 0; j < domainPartition + 1; j++)    // Y-Axis
    {
      for (auto i = 0; i < domainPartition + 1; i++)  // X-Axis
      {
        feature.emplace_back( (
          functionSet.*fC.indexedFunc[initFuncInd])(i*xStep, j*yStep, k*zStep) );
      }
    }
  }
}
//=================================================================================================