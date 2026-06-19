#pragma once

#include <math.h>

constexpr int D = 2;
constexpr int Q = 9;

extern __constant__ float w[Q];
extern __constant__ float c_ix[Q];
extern __constant__ float c_iy[Q];

const float a_s = sqrtf(3);
const float a_s2 = 3.f;
const float a_s4 = 9.f;
const float inv_as2 = 1.f / a_s2;