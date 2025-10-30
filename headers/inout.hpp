#include "settings.hpp"
/*! \brief input function
 *  \param[in] filename - the name of input file
 *  \param[in] params   - the link to parameters structure
 *  \return code of operation
 */
int input(const std::string& filename, model_data& params);

/*! \brief mesh output function
 *  \param[in] filename - the name of input file
 *  \param[in] axisMesh - mesh along <AXIS>
 *  \return code of operation
 */
int meshOutput(const std::string& filename, const std::vector<double>& axisMesh);