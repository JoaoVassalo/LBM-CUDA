#include "stencil.cuh"

extern __constant__ float w[Q] = {
    4.f/9.f,
    1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f,
    1.f/36.f, 1.f/36.f, 1.f/36.f, 1.f/36.f
};
extern __constant__ int c_ix[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
extern __constant__ int c_iy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};