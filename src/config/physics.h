#pragma once

#include "geometry.h"
#include "stencil.cuh"

// Calc Tau.
constexpr int Re = 1000;
constexpr float u_max = 0.0256f;
constexpr int delta_t = 1;
constexpr float ni = u_max * (float)Ny / (float)Re;
const float tau = ni * a_s2 + 0.5f;
const float omega = 1.0f / tau;