#pragma once

#include "../../config/geometry.h"

#include "build_mom.cuh"

#include <cstdint>

struct Grid2D
{
    uint8_t *mask;
    uint8_t *node;

    static constexpr int size = Nx * Ny * sizeof(uint8_t);
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
        break;
    case 1:
        return 3;
        break;
    case 2:
        return 4;
        break;
    case 3:
        return 1;
        break;
    case 4:
        return 2;
        break;
    case 5:
        return 7;
        break;
    case 6:
        return 8;
        break;
    case 7:
        return 5;
        break;
    default: // case 8
        return 6;
        break;
    }
}

__host__ void build_grid(D2Q9 sim, Grid2D &grid);

__global__ void init_grid(uint8_t *node, uint8_t *mask);