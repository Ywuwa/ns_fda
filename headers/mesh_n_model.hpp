#include "settings.hpp"
#include <functional>
//=================================== INITIAL FUCNTIONS ===========================================
/*! \brief Function set consits of 3 functions for the velocity components and the pressure
 */
struct functionSet
{
  double initV1(const double x, const double y, const double z);
  double initV2(const double x, const double y, const double z);
  double initV3(const double x, const double y, const double z);
  double initPress(const double x, const double y, const double z);  
};
/*! \brief Function Container accumulate 4 functions of chosen function Set
 */
template <typename fSet>
struct functionContainer
{
  std::vector<double (*)(const double, const double, const double)> indexedFunc;
  functionContainer(fSet& name)
  {
    indexedFunc.emplace_back(fSet::initV1);
    indexedFunc.emplace_back(fSet::initV2);
    indexedFunc.emplace_back(fSet::initV3);
    indexedFunc.emplace_back(fSet::initPress);
  }
};
//=================================================================================================

//=================================== INITIAL CONDITIONS ==========================================
/*! \brief Applying initial conditions for functions
 *  \param[in] feature     - feature that needs to be initialized
 *  \param[in] initFuncInd - the required function index (0, 1, 2 or 3)
 *  \param[in] params      - the link to parameters structure
 */
void initialConditions(
  std::vector<double>& feature, const uint initFuncInd, const model_data& params);
//=================================================================================================

//===================================== AXIS MESH INIT ============================================
/*! \brief axis mesh init function
 *  \param[in] axisMesh - mesh along <AXIS>
 *  \param[in] params   - the link to parameters structure
 *  \param[in] axisLen  - the length of the current axis
 */
void meshInit(std::vector<double>& axisMesh, const model_data& params, const double axisLen);
//=================================================================================================