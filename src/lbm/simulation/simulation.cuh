#pragma once

#include "../../config/gridConfig.cuh"
#include "../../config/momConfig.cuh"
#include "../../config/simulationConfig.cuh"
#include "../../config/outputConfig.cuh"

#include "../../io/writeOutput.cuh"

#include "layerPipeline.cuh"

#include <chrono>

__host__ void simulation(moments sim, Grid2D grid);