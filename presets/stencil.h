#pragma once

#include <math.h>

constexpr int D = 2;
constexpr int Q = 9;
__constant__ float w[Q] = {
    4.f/9.f,
    1.f/9.f, 1.f/9.f, 1.f/9.f, 1.f/9.f,
    1.f/36.f, 1.f/36.f, 1.f/36.f, 1.f/36.f
};
__constant__ int c_ix[Q] = {0, 1, 0, -1, 0, 1, -1, -1, 1};
__constant__ int c_iy[Q] = {0, 0, 1, 0, -1, 1, 1, -1, -1};
const float a_s = sqrtf(3);
const float a_s2 = 3.f;
const float a_s4 = 9.f;
const float inv_as2 = 1.f/a_s2;