#pragma once

#include "../presets/config.h"

__global__ void lbm_step(varUnit *mom_in, varUnit *mom_out);