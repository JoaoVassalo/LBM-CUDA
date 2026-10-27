#pragma once

#define varType float

#include "stencilConfig.cuh"

namespace timeConfig
{
    inline constexpr size_t tf = 10000;
    inline constexpr size_t printNumber = 10;
    inline constexpr size_t tInterval = tf / printNumber;
}

namespace physics
{
    inline constexpr int Re = 1000;
    inline constexpr varType u_max = 0.0256;
    inline constexpr int delta_t = 1;
    inline constexpr varType ni = u_max * (varType)Geometry::NY / (varType)Re;
    inline constexpr varType tau = ni * D2Q9::a_s2 + 0.5;
    inline constexpr varType omega = 1.0 / tau;
}
