#include <fstream>
#include <string>

#include "settings.hpp"
/*! \brief input function
 *  \param[in] filename - the name of input file
 *  \param[in] params   - the link to parameters structure
 *  \return code of operation
 */
int dataInput(const std::string& filename, model_data& params);

/*! \brief mesh output function
 *  \param[in] filename - the name of input file
 *  \param[in] axisMesh - mesh along <AXIS>
 *  \return code of operation
 */
int meshOutput(const std::string& filename, const std::vector<double>& axisMesh);

/*! \brief grid-function value output
 *  \param[in] path      - path to dir
 *  \param[in] name      - direct filename
 *  \param[in] extension - file extension
 *  \param[in] funcValue - function value
 *  \param[in] params    - the link to parameters structure
 *  \return code of operation
 */
int funcOutput(const std::string& path, const std::string& name, const std::string& extension, 
               const std::vector<double>& funcValue, const model_data& params);
  //=================================================================================================