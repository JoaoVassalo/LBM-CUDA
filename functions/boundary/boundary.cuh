#pragma once

#include "../../cuda_config/var.cuh"

#include "../../presets/config.h"
#include "../../presets/physics.h"
#include "../../presets/geometry.h"
#include "../../presets/stencil.cuh"

#include "../equations/f_i.cuh"

#include "../grid_id.cuh"

__device__ void boundary(int size, CInt *Is, CInt *Os, int x, int y, varUnit *mom_in,
                         varUnit *mom_out);