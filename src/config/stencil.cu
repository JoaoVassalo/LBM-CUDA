#include "stencil.cuh"

__constant__ float w[Q] = {
    4.f / 9.f,
    1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f,
    1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f};
__constant__ float c_ix[Q] = {0.f, 1.f, 0.f, -1.f, 0.f, 1.f, -1.f, -1.f, 1.f};
__constant__ float c_iy[Q] = {0.f, 0.f, 1.f, 0.f, -1.f, 1.f, 1.f, -1.f, -1.f};