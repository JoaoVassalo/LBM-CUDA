#pragma once

#include "../../presets/config.h"
#include "../../presets/stencil.cuh"
#include "../../presets/physics.h"
#include "../../cuda_config/var.cuh"
#include "../from_id.cuh"
#include "../grid_id.cuh"
#include "../equations/f_i.cuh"

__device__ void inlet_north(CInt *I_s, int x, int y,
                            float *mom_in,
                            float *mom_out);

__device__ void inlet_south(CInt *I_s, int x, int y,
                            float *mom_in,
                            float *mom_out);

__device__ void inlet_east(CInt *I_s, int x, int y,
                           float *mom_in,
                           float *mom_out);

__device__ void inlet_west(CInt *I_s, int x, int y,
                           float *mom_in,
                           float *mom_out);

__device__ void inlet_northeast(CInt *I_s, int x, int y,
                                float *mom_in,
                                float *mom_out);

__device__ void inlet_northwest(CInt *I_s, int x, int y,
                                float *mom_in,
                                float *mom_out);

__device__ void inlet_southeast(CInt *I_s, int x, int y,
                                float *mom_in,
                                float *mom_out);

__device__ void inlet_southwest(CInt *I_s, int x, int y,
                                float *mom_in,
                                float *mom_out);