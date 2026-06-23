#pragma once

#include "vtk.cuh"

#include "../config/momConfig.cuh"
#include "../config/gridConfig.cuh"

#include <string>

__host__ void writeOutput(moments sim, Grid2D grid, size_t t, size_t momByteSize, const std::string path);