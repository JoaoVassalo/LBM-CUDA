#pragma once

#include "../../presets/config.h"

__device__ void wall(CInt *I_s, CInt *O_s, int x, int y,
                     float *mom_in,
                     float *mom_out);