#include "D:\Projects\C++\ns_fda\headers\compute_flow.hpp"

//==================================== FLOW COMPUTATION ===========================================
void compute(
  model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p)
{
  uint tick = 0;
  const auto dimSize ( params.domainPartition );
  const double timeStep ( params.duration / params.timePartition );
  const double hX (params.xLen / dimSize );
  const double hY (params.yLen / dimSize );
  const double hZ (params.zLen / dimSize );

  const size_t vecSize = (dimSize + 1) * (dimSize + 1) * (dimSize + 1);
  std::vector<double> u1(vecSize);
  std::vector<double> v1(vecSize);
  std::vector<double> w1(vecSize);
  
  while (tick < params.timePartition + 1)
  {
    // velocity compute
    //---------------------------- inner knots --------------------------------
    for (auto k = 1; k < dimSize; k++)      // Z-Axis
    {
      for (auto j = 1; j < dimSize; j++)    // Y-Axis
      {
        for (auto i = 1; i < dimSize; i++)  // X-Axis
        {
          uint index (k*(dimSize+1)*(dimSize+1) + j*(dimSize+1) + i);
          u1[index] = u[index] - timeStep;
          v1[index] = v[index] - timeStep;
          w1[index] = w[index] - timeStep;
        }
      }
    }
    //-------------------------------------------------------------------------
    //---------------------------- border knots -------------------------------
    // XY-plane Z = 0 / Z = MAX
    for (auto j = 1; j < dimSize; j++)    // Y-Axis
    {
      for (auto i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index1 (j*(dimSize+1) + i);
        uint index2 (dimSize*(dimSize+1)*(dimSize+1) + j*(dimSize+1) + i);
      }
    }
    // XZ-plane Y = 0 / Y = MAX
    for (auto k = 1; k < dimSize; k++)    // Z-Axis
    {
      for (auto i = 1; i < dimSize; i++)  // X-Axis
      {
        uint index1 (k*(dimSize+1)*(dimSize+1) + i);
        uint index2 (k*(dimSize+1)*(dimSize+1) + dimSize*(dimSize+1) + i);
      }
    }
    // YZ-plane X = 0 / X = MAX
    for (auto k = 1; k < dimSize; k++)    // Z-Axis
    {
      for (auto j = 1; j < dimSize; j++)  // X-Axis
      {
        uint index1 (k*(dimSize+1)*(dimSize+1) + j*(dimSize+1));
        uint index2 (k*(dimSize+1)*(dimSize+1) + j*(dimSize+1) + dimSize);
      }
    }
    //-------------------------------------------------------------------------
    //---------------------------- edge knots ---------------------------------
    // Y = Z = (0 || MAX)
    for (auto i = 1; i < dimSize; i++)
    {
      // Y = Z = 0
      uint index1 (i);
      // Y = MAX, Z = 0
      uint index2 (dimSize*(dimSize+1) + i);
      // Y = 0, Z = MAX
      uint index3 (dimSize*(dimSize+1)*(dimSize+1) + i);
      // Y = MAX, Z = MAX
      uint index4 (dimSize*(dimSize+1)*(dimSize+1) + dimSize*(dimSize+1) + i);
    }
    // X = Z = (0 || MAX)
    for (auto j = 1; j < dimSize; j++)
    {
      // X = Z = 0
      uint index1 (j*(dimSize+1));
      // X = MAX, Z = 0
      uint index2 (j*(dimSize+1) + dimSize);
      // X = 0, Z = MAX
      uint index3 (dimSize*(dimSize+1)*(dimSize+1) + j*(dimSize+1));
      // X = MAX, Z = MAX
      uint index4 (dimSize*(dimSize+1)*(dimSize+1) + j*(dimSize+1) + dimSize);
    }
    // X = Y = (0 || MAX)
    for (auto k = 1; k < dimSize; k++)
    {
      // X = Y = 0
      uint index1 (k*(dimSize+1)*(dimSize+1));
      // X = MAX, Y = 0
      uint index2 (k*(dimSize+1)*(dimSize+1) + dimSize);
      // X = 0, Y = MAX
      uint index3 (k*(dimSize+1)*(dimSize+1) + dimSize*(dimSize+1));
      // X = MAX, Y = MAX
      uint index4 (k*(dimSize+1)*(dimSize+1) + dimSize*(dimSize+1) + dimSize);
    }
    //-------------------------------------------------------------------------
    //------------------------- cube vertices ---------------------------------
    // X = Y = Z = 0
    uint index (0);
    // X = MAX, Y = Z = 0
    index = dimSize;
    // X = Z = 0, Y = MAX
    index = dimSize*(dimSize+1);
    // X = Y = 0, Z = MAX
    index = dimSize*(dimSize+1)*(dimSize+1);
    // X = Y = MAX, Z = 0
    index = dimSize*(dimSize+1) + dimSize;
    // X = Z = MAX, Y = 0
    index = dimSize*(dimSize+1)*(dimSize+1) + dimSize;
    // Y = Z = MAX, X = 0
    index = dimSize*(dimSize+1)*(dimSize+1) + dimSize*(dimSize+1);
    // X = Y = Z = MAX
    index = dimSize*(dimSize+1)*(dimSize+1) + dimSize*(dimSize+1) + dimSize;
    //-------------------------------------------------------------------------
    // broadcast, sync point

    // pressure compute
    tick += 1;
  }
  std::cout << tick << '\n';
}
//=================================================================================================