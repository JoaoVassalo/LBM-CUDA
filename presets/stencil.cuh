#pragma once

#include <math.h>

constexpr int D = 2;
constexpr int Q = 9;

extern __device__ __constant__ float w[Q];
extern __device__ __constant__ float c_ix[Q];
extern __device__ __constant__ float c_iy[Q];

extern __constant__ int Is_SW[4];
extern __constant__ int Os_SW[4];

extern __constant__ int Is_SE[4];
extern __constant__ int Os_SE[4];

extern __constant__ int Is_NW[4];
extern __constant__ int Os_NW[4];

extern __constant__ int Is_NE[4];
extern __constant__ int Os_NE[4];

extern __constant__ int Is_N[6];
extern __constant__ int Os_N[6];

extern __constant__ int Is_S[6];
extern __constant__ int Os_S[6];

extern __constant__ int Is_W[6];
extern __constant__ int Os_W[6];

extern __constant__ int Is_E[6];
extern __constant__ int Os_E[6];

const float a_s = sqrtf(3);
const float a_s2 = 3.f;
const float a_s4 = 9.f;
const float inv_as2 = 1.f / a_s2;