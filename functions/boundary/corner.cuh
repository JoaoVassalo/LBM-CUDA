#pragma once

#include "../../presets/config.h"

__device__ void corner(CInt *I_s, CInt *O_s, int x, int y, float *rho_in, float *ux_in, float *uy_in, float *mxx_in, float *mxy_in, float *myy_in,
                       float *rho_out);