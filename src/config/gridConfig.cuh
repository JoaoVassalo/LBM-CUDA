#pragma once

#include "stencilConfig.cuh"

#include <cstdint>

struct Grid2D
{
    uint8_t *mask;
    uint8_t *node;

    static constexpr int gridByteSize = Geometry::NX * Geometry::NY * sizeof(uint8_t);
};
