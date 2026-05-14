#pragma once

#include <cstdint>

__device__ void applyBoundary(uint8_t *node, uint8_t *mask, float *mom_in, float *mom_out);
