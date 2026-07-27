#pragma once

#include "../../presets/config.h"
#include "../../presets/stencil.cuh"
#include "../../cuda_config/var.cuh"
#include "../from_id.cuh"
#include "../grid_id.cuh"
#include "../equations/f_i.cuh"

__device__ void outlet_north(CInt *I_s, int x, int y,
                             float *mom_in,
                             float *mom_out);

__device__ void outlet_south(CInt *I_s, int x, int y,
                             float *mom_in,
                             float *mom_out);

__device__ void outlet_east(CInt *I_s, int x, int y,
                            float *mom_in,
                            float *mom_out);

__device__ void outlet_west(CInt *I_s, int x, int y,
                            float *mom_in,
                            float *mom_out);

__device__ void outlet_northeast(CInt *I_s, int x, int y,
                                 float *mom_in,
                                 float *mom_out);

__device__ void outlet_northwest(CInt *I_s, int x, int y,
                                 float *mom_in,
                                 float *mom_out);

__device__ void outlet_southeast(CInt *I_s, int x, int y,
                                 float *mom_in,
                                 float *mom_out);

__device__ void outlet_southwest(CInt *I_s, int x, int y,
                                 float *mom_in,
                                 float *mom_out);