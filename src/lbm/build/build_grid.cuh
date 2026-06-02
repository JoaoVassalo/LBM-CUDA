#pragma once

#include <cstdint>
#include "../../config/geometry.h"
#include "build_mom.cuh"

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

__host__ __device__ __forceinline__ int income(int i)
{
    switch (i)
    {
    case 0:
        return 0;
    case 1:
        return 3;
    case 2:
        return 4;
    case 3:
        return 1;
    case 4:
        return 2;
    case 5:
        return 7;
    case 6:
        return 8;
    case 7:
        return 5;
    default:
        return 6;
    }
}

__host__ void build_grid(D2Q9 sim, Grid2D &grid);

__global__ void init_grid(uint8_t *node, uint8_t *mask);
