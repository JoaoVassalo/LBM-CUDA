#pragma once

#include "../config/momConfig.cuh"
#include "../config/stencilConfig.cuh"
#include "../core/indexing.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>

__host__ void writeVTI(size_t t, const std::string &outDir, float *mom_host);