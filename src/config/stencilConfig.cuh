#pragma once

#include <cstddef>
#include <math.h>

namespace Geometry
{
    // Geometry definition
    constexpr int NX = 2048;
    constexpr int NY = 2048;
    constexpr int NZ = 0;

    // Layer definition
    constexpr int LNX = NX;
    constexpr int LNY = 3;
    constexpr int LNZ = 0;
};

// Defining block size and number for kernel initialization, in this branch we use this part mainly for the initialization part.
#define BX 32
#define BY 16
#define GX ((Geometry::NX + BX - 1) / BX)
#define GY ((Geometry::NY + BY - 1) / BY)

inline dim3 block(BX, BY);
inline dim3 blockNumber(GX, GY);

// Defining block size and number for kernel using layer.
constexpr int layerBlockX = Geometry::LNX < 256 ? Geometry::LNX : 256;

inline dim3 blockLayer = dim3(layerBlockX, 1);
inline dim3 blockNumberLayer = dim3((Geometry::LNX + layerBlockX - 1) / layerBlockX, 1);

// Stencil definition
namespace D2Q9
{
    constexpr int D = 2;
    constexpr int Q = 9;

    constexpr int momNum = 6; // Moment number for D2Q9 is set at 6. Look at momConfig.cuh for more details.

    constexpr size_t momSize = static_cast<size_t>(Geometry::NX) * Geometry::NY;
    constexpr size_t momByteSize = momSize * momNum * sizeof(varType);
    constexpr size_t outputByteSize = momSize * 6 * sizeof(varType);

    constexpr size_t layerSize = Geometry::NX;
    constexpr size_t layerByteSize = layerSize * momNum * sizeof(varType);

    const varType a_s = sqrt(static_cast<varType>(3));
    constexpr varType a_s2 = 3.;
    constexpr varType a_s4 = 9.;
    constexpr varType inv_as2 = 1. / a_s2;

    __host__ __device__ __forceinline__ varType w(int i)
    {
        switch (i)
        {
        case 0:
            return 4. / 9.;
        case 1:
        case 2:
        case 3:
        case 4:
            return 1. / 9.;
        case 5:
        case 6:
        case 7:
        case 8:
            return 1. / 36.;
        default:
            return 0.;
        }
    }

    __host__ __device__ __forceinline__ varType cx(int i)
    {
        switch (i)
        {
        case 0:
            return 0.;
        case 1:
            return 1.;
        case 2:
            return 0.;
        case 3:
            return -1.;
        case 4:
            return 0.;
        case 5:
            return 1.;
        case 6:
            return -1.;
        case 7:
            return -1.;
        case 8:
            return 1.;
        default:
            return 0.;
        }
    }

    __host__ __device__ __forceinline__ varType cy(int i)
    {
        switch (i)
        {
        case 0:
            return 0.;
        case 1:
            return 0.;
        case 2:
            return 1.;
        case 3:
            return 0.;
        case 4:
            return -1.;
        case 5:
            return 1.;
        case 6:
            return 1.;
        case 7:
            return -1.;
        case 8:
            return -1.;
        default:
            return 0.;
        }
    }

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
};
