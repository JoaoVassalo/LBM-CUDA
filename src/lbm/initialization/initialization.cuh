#pragma once

#include "initDomain.cuh"
#include "initGrid.cuh"

#include "../../config/momConfig.cuh"
#include "../../config/gridConfig.cuh"
#include "../../config/outputConfig.cuh"

#include "../../io/writeOutput.cuh"

#include <iostream>
#include <string>
#include <fstream>
#include <cstdint>
#include <bitset>

__host__ void initialization(moments &sim, Grid2D &grid);