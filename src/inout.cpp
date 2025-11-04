#include "D:\Projects\C++\ns_fda\headers\inout.hpp"

//====================================== FILE READING =============================================
int dataInput(const std::string& filename, model_data& params) {
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
	std::vector<std::string> data_collection;
	while (std::getline(file, line)) {
		if (!(line.empty())) {
			size_t colon_pos = line.find(':');
			if (colon_pos == std::string::npos || colon_pos == line.length() - 1) {
				continue;  // Skip lines that don't contain "text:number"
			}

			std::string number_str = line.substr(colon_pos + 1);
			data_collection.emplace_back(number_str.begin(), number_str.end());
		}
	}
	// We parse this file under the assumption, that inputFile configuration is known
	// If new parameters are added into inputFile, the code below have to be modified
	params.domainPartition = std::stoi(data_collection[0]);
	params.xLen = std::stof(data_collection[1]);
	params.yLen = std::stof(data_collection[2]);
	params.zLen = std::stof(data_collection[3]);
	params.Reyn = std::stof(data_collection[4]);
	return 0;
}
//=================================================================================================

//======================================= MESH OUTPUT =============================================
int meshOutput(const std::string& filename, const std::vector<double>& axisMesh) {
	// Open a file for writing, std::ios::out ensures it's opened for writing.
	// If the file doesn't exist, it will be created. If it exists, its content will be overwritten
	std::ofstream outputFile(filename, std::ios::out);

	// Check if the file was opened successfully
	if (outputFile.is_open()) {
		// Iterate through the array and write each element to the file
		for (auto knot : axisMesh) {
				outputFile << knot << std::endl; // Write the number and a newline character
		}
		// Close the file stream
		outputFile.close();
		std::cout << "mesh successfully written to " << filename << std::endl;
	} else {
		std::cerr << "Error: Unable to open file for writing." << std::endl;
	}

	return 0;
}
//=================================================================================================

//===================================== FUNCTION OUTPUT ===========================================
int funcOutput(const std::string& path, const std::string& name, const std::string& extension, 
               const std::vector<double>& funcValue, const model_data& params)
{
	const std::string filename = path + name + extension;
  std::ofstream outputFile(filename, std::ios::out);
	size_t index  (0);
	// Check if the file was opened successfully
	if (outputFile.is_open()) {
		const uint domainPartition = params.domainPartition;
		for (auto k = 0; k < domainPartition + 1; k++)      // Z-Axis
    {
			for (auto j = 0; j < domainPartition + 1; j++)    // Y-Axis
			{
				for (auto i = 0; i < domainPartition + 1; i++)  // X-Axis
				{
					index = k*(domainPartition+1)*(domainPartition+1) + j*(domainPartition+1) + i;
					outputFile << funcValue[index] << ' ';
				}
				outputFile << '\n';
			}
			outputFile << '\n';
		}
		// Close the file stream
		outputFile.close();
		std::cout << name << " successfully written to " << filename << std::endl;
	} else {
		std::cerr << "Error: Unable to open file for writing." << std::endl;
	}

	return 0;
}
//=================================================================================================