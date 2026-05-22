#pragma once

#include <cstdint>

__device__ void center(int x, int y,
                       D2Q9 sim);

__device__ void boundary(uint8_t mask_in, uint8_t node_in, int x, int y,
                         D2Q9 sim);
