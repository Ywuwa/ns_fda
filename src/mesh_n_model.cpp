#include "D:\Projects\C++\ns_fda\headers\mesh_n_model.hpp"

//=================================== INITIAL FUCNTIONS ===========================================
double functionSet::initV1(const double x, const double y, const double z)
{
  return std::max(x, y);
}
double functionSet::initV2(const double x, const double y, const double z)
{
  return 1.0;
}
double functionSet::initV3(const double x, const double y, const double z)
{
  return 2.0;
}
double functionSet::initPress(const double x, const double y, const double z)
{
  return 3.0;
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
  functionSet fSet;
  functionContainer fC(fSet);

  for (auto k = 0; k < domainPartition + 1; k++)      // Z-Axis
  {
    for (auto j = 0; j < domainPartition + 1; j++)    // Y-Axis
    {
      for (auto i = 0; i < domainPartition + 1; i++)  // X-Axis
      {
        feature.emplace_back( (fSet.*fC.indexedFunc[initFuncInd])(i*xStep, j*yStep, k*zStep) );
      }
    }
  }
}
//=================================================================================================