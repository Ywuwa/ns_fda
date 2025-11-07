#include "settings.hpp"

//==================================== FLOW COMPUTATION ===========================================
/*! \brief Flow computation function
 *  \param[in] params - model data
 *  \param[in] u, v, w - velocity components
 *  \param[in] p - pressure
 */
void compute(
  model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p);
//=================================================================================================