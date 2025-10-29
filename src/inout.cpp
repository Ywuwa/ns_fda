#include "D:\Projects\C++\ns_fda\headers\inout.hpp"
//====================================== FILE READING =============================================

int input(const std::string& filename, model_data& params) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        return -1; // Cannot open file
    }

    // Check if file is empty
    file.peek();
    if (file.eof()) {
        return -2; // File is empty
    }

    std::string line;
    // Skip first 2 intro lines
    for (int i = 0; i < 2; ++i) {
        if (!std::getline(file, line)) {
            return -3; // Cannot read 3 lines
        }
    }
    std::cout << "before unboxing" << '\n';
    std::vector<std::string> data_collection;
    while (std::getline(file, line)) {
      std::cout << line << '\n';
      if (!(line.empty())) {
          size_t colon_pos = line.find(':');
          if (colon_pos == std::string::npos || colon_pos == line.length() - 1) {
              continue;  // Skip lines that don't contain "text:number"
          }

          std::string number_str = line.substr(colon_pos + 1);
          std::cout << number_str << '\n';
          data_collection.emplace_back(number_str.begin(), number_str.end());
      }
      std::cout << line << '\n';
    }
    // We parse this file under the assumption, that inputFile configuration is known
    // If new parameters are added into inputFile, the code below have to be modified
    params.mesh_knots_num = std::stoi(data_collection[0]);
    params.Reyn = std::stof(data_collection[1]);
    return 0;
}
//=================================================================================================