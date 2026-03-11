#pragma once

#include "/home/jvassalo/LBM-CUDA/presets/geometry.h"
#include "/home/jvassalo/LBM-CUDA/presets/stencil.cuh"


//Calc Tau.
constexpr int Re = 1000;
constexpr float u_max = 0.0256;
constexpr int delta_t = 1;
constexpr float ni = u_max*Ny/Re;
const float tau = ni*a_s2;
const float omega = 1/tau;

