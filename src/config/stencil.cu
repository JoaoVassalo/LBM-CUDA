#include "stencil.cuh"

extern __constant__ float w[Q] = {
    4.f / 9.f,
    1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f,
    1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f};
extern __constant__ float c_ix[Q] = {0.f, 1.f, 0.f, -1.f, 0.f, 1.f, -1.f, -1.f, 1.f};
extern __constant__ float c_iy[Q] = {0.f, 0.f, 1.f, 0.f, -1.f, 1.f, 1.f, -1.f, -1.f};
