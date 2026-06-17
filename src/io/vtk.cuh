#pragma once

#include "../config/geometry.h"
#include "../core/indexing.cuh"

#include "../lbm/build/build_mom.cuh"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <iostream>
#include <string>

__host__ void write_vti(int step, const std::string &out_dir, float *mom_host);
