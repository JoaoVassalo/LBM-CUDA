#pragma once

#include "../../presets/physics.h"
#include "../../presets/geometry.h"

__device__ __forceinline__ varUnit velocity(int y)
{
    return (((-4.f * u_max) / (Ny * Ny)) * ((varUnit)y * (varUnit)y) + ((4.f * u_max) / (Ny)) * (varUnit)y);
}