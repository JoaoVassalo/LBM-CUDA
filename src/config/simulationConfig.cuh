#pragma once

#include "stencilConfig.cuh"

namespace time
{
    size_t tf = 4e6;
    size_t printNumber = 100;
    size_t tInterval = tf / printNumber;
}

namespace physics
{
    constexpr int Re = 1000;
    constexpr float u_max = 0.0256f;
    constexpr int delta_t = 1;
    constexpr float ni = u_max * (float)Geometry::NY / (float)Re;
    const float tau = ni * a_s2 + 0.5f;
    const float omega = 1.0f / tau;
}