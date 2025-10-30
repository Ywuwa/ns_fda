#include "D:\Projects\C++\ns_fda\headers\mesh_n_model.hpp"

//===================================== AXIS MESH INIT ============================================
void meshInit(std::vector<double>& axisMesh, const model_data& params, const double axisLen)
{
  const double axisStep = axisLen / params.domainPartition;
  axisMesh.reserve(params.domainPartition + 1);
  for (size_t i = 0; i < params.domainPartition + 1; i++)
  {
    axisMesh[i] = i * axisStep;
  }
}
//=================================================================================================