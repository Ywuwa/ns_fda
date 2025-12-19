#include "settings.hpp"

//============================== FLOW COMPUTATION via FDA1/FDA3 ===================================
/*! \brief Flow computation function in the cubic domain, based on FDA1 or FDA3
 *  \param[in] params - model data
 *  \param[in] u, v, w - velocity components
 *  \param[in] p - pressure
 */
void compute_cube_FDA1_3(
  const model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p);
//=================================================================================================

//============================== FLOW COMPUTATION via FDA2/FDA4 ===================================
/*! \brief Flow computation function in the cubic domain, based on FDA2 or FDA4
 *  \param[in] params - model data
 *  \param[in] u, v, w - velocity components
 *  \param[in] p - pressure
 */
void compute_cube_FDA2_4(
  const model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p0);
//=================================================================================================

//======================================== RESIDUALS ==============================================
/*! \brief Velocity residual computation function, based on FDA1 and FDA3
 *  \param[in] params - model data
 *  \param[in] uEst, vEst, wEst - estimated velocity components
 *  \param[in] pEst - estimated pressure
 *  \param[in] uExac, vExac, wExac - exact velocity components
 *  \param[in] pExacc - exact pressure
 *  \return residial value
 */
double velocity_residual_FDA1_3(
  const model_data& params, 
  std::vector<double>& uEst, std::vector<double>& vEst, std::vector<double>& wEst, 
  std::vector<double>& pEst,
  std::vector<double>& uExac, std::vector<double>& vExac, std::vector<double>& wExac, 
  std::vector<double>& pExac);

  /*! \brief Velocity residual computation function, based on FDA2 and FDA4
 *  \param[in] params - model data
 *  \param[in] uEst, vEst, wEst - estimated velocity components
 *  \param[in] pEst - estimated pressure
 *  \param[in] uExac, vExac, wExac - exact velocity components
 *  \param[in] pExacc - exact pressure
 *  \return residial value
 */
double velocity_residual_FDA2_4(
  const model_data& params, 
  std::vector<double>& uEst, std::vector<double>& vEst, std::vector<double>& wEst, 
  std::vector<double>& pEst,
  std::vector<double>& uExac, std::vector<double>& vExac, std::vector<double>& wExac, 
  std::vector<double>& pExac);
//=================================================================================================