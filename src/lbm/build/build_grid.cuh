#pragma once

#include <cstdint>
#include "../../config/geometry.h"

struct Grid2D
{
    uint8_t *mask;
    uint8_t *node;

    static constexpr int size = Geometry::Nx * Geometry::Ny * sizeof(uint8_t);
};

enum class domainTags : uint8_t
{
    Fluid = 0,
    Boundary = 1,
};

__host__ __device_builtin__ __forceinline__ int income(int i);

__host__ void build_grid(D2Q9 sim, Grid2D grid);

__global__ void init_grid(uint8_t *node, uint8_t *mask);