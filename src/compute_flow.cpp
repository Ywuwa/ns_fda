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
  while (tick < params.timePartition + 1)
  {
    // velocity compute

    // broadcast, sync point

    // pressure compute
    tick += 1;
  }
  std::cout << tick << '\n';
}
//=================================================================================================