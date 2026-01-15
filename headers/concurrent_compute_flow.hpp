#include "settings.hpp"
#include <mutex>
#include <thread>

struct ConcurrentSharedData {
    std::mutex mtx;
    model_data params;

    std::vector<double> u;
    std::vector<double> v;
    std::vector<double> w;
    std::vector<double> p;

    std::vector<double> uHat;
    std::vector<double> vHat;
    std::vector<double> wHat;
    std::vector<double> pHat;

    std::vector<double> uExac;
    std::vector<double> vExac;
    std::vector<double> wExac;
    std::vector<double> pExac;
    ConcurrentSharedData(
      model_data&& model, std::vector<double>&& pIn,
      std::vector<double>&& uIn, std::vector<double>&& vIn, std::vector<double>&& wIn
    ) 
    {
      params = std::move(model);
      u = std::move(uIn); v = std::move(vIn); u = std::move(vIn); p = std::move(pIn);
      uHat.resize(u.size()); vHat.resize(v.size()); wHat.resize(w.size()); pHat.resize(p.size());
      uExac.resize(u.size()); vExac.resize(v.size()); wExac.resize(w.size()); pExac.resize(p.size());
    }

};

//============================== FLOW COMPUTATION via FDA1/FDA3 ===================================
/*! \brief Flow computation function in the cubic domain, based on FDA1 or FDA3, CONCURRENT
 *  \param[in] params - model data
 *  \param[in] u, v, w - velocity components
 *  \param[in] p - pressure
 */
void conc_compute_cube_FDA1_3( const model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p);
//=================================================================================================

//============================== FLOW COMPUTATION via FDA2/FDA4 ===================================
/*! \brief Flow computation function in the cubic domain, based on FDA2 or FDA4, CONCURRENT
 *  \param[in] params - model data
 *  \param[in] u, v, w - velocity components
 *  \param[in] p - pressure
 */
void conc_compute_cube_FDA2_4( const model_data& params, 
  std::vector<double>& u, std::vector<double>& v, std::vector<double>& w, std::vector<double>& p);
//=================================================================================================

//======================================== RESIDUALS ==============================================
/*! \brief Velocity residual computation function, based on FDA1 and FDA3, CONCURRENT
 *  \param[in] params - model data
 *  \param[in] uEst, vEst, wEst - estimated velocity components
 *  \param[in] pEst - estimated pressure
 *  \param[in] uExac, vExac, wExac - exact velocity components
 *  \param[in] pExacc - exact pressure
 *  \return residial value
 */
double conc_velocity_residual_FDA1_3(
  const model_data& params, 
  std::vector<double>& uEst, std::vector<double>& vEst, std::vector<double>& wEst, 
  std::vector<double>& pEst,
  std::vector<double>& uExac, std::vector<double>& vExac, std::vector<double>& wExac, 
  std::vector<double>& pExac);

  /*! \brief Velocity residual computation function, based on FDA2 and FDA4, CONCURRENT
 *  \param[in] params - model data
 *  \param[in] uEst, vEst, wEst - estimated velocity components
 *  \param[in] pEst - estimated pressure
 *  \param[in] uExac, vExac, wExac - exact velocity components
 *  \param[in] pExacc - exact pressure
 *  \return residial value
 */
double conc_velocity_residual_FDA2_4(
  const model_data& params, 
  std::vector<double>& uEst, std::vector<double>& vEst, std::vector<double>& wEst, 
  std::vector<double>& pEst,
  std::vector<double>& uExac, std::vector<double>& vExac, std::vector<double>& wExac, 
  std::vector<double>& pExac);
//=================================================================================================