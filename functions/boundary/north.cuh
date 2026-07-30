#pragma once

#include "../../presets/config.h"

__device__ void north(int size, CInt *I_s, CInt *O_s, int x, int y,
                      varUnit *mom_in,
                      varUnit *mom_out);