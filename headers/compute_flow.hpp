#include "settings.hpp"

//==================================== FLOW COMPUTATION ===========================================
/*! \brief Flow computation function
 *  \param[in] params - model data
 *  \param[in] u, v, w - velocity components
 *  \param[in] p - pressure
 */
void compute(
  const model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p);
//=================================================================================================

//======================================== RESIDUALS ==============================================
/*! \brief Flow computation function
 *  \param[in] params - model data
 *  \param[in] uEst, vEst, wEst - estimated velocity components
 *  \param[in] pEst - estimated pressure
 *  \param[in] uExac, vExac, wExac - exact velocity components
 *  \param[in] pExacc - exact pressure
 *  \return residial value
 */
double velocity_residual(
  const model_data& params, 
  std::vector<double>& uEst, std::vector<double>& vEst, std::vector<double>& wEst, 
  std::vector<double>& pEst,
  std::vector<double>& uExac, std::vector<double>& vExac, std::vector<double>& wExac, 
  std::vector<double>& pExac);
//=================================================================================================