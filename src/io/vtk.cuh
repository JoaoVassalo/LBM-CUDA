#pragma once

#include "../config/momConfig.cuh"
#include "../config/simulationConfig.cuh"
#include "../core/indexing.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>

__host__ void writeVTI(size_t t, const std::string &outDir, varType *mom_host);
