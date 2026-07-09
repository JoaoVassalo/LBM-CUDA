#pragma once

#include "stencilConfig.cuh"

namespace timeConfig
{
    inline constexpr size_t tf = 100;
    inline constexpr size_t printNumber = 100;
    inline constexpr size_t tInterval = tf / printNumber;
}

namespace physics
{
    inline constexpr int Re = 1000;
    inline constexpr float u_max = 0.0256f;
    inline constexpr int delta_t = 1;
    inline constexpr float ni = u_max * (float)Geometry::NY / (float)Re;
    inline constexpr float tau = ni * D2Q9::a_s2 + 0.5f;
    inline constexpr float omega = 1.0f / tau;
}