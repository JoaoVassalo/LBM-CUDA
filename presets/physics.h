#pragma once

#include "/home/jvassalo/LBM-CUDA/presets/geometry.h"
#include "/home/jvassalo/LBM-CUDA/presets/stencil.cuh"

// Calc Tau.
constexpr int Re = 100;
constexpr varUnit u_max = 0.0256f;
constexpr int delta_t = 1;
constexpr varUnit ni = u_max * (varUnit)Ny / (varUnit)Re;
const varUnit tau = ni * a_s2 + 0.5f;
const varUnit omega = 1.0f / tau;