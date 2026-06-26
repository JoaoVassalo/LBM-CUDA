#pragma once

#include <math.h>

namespace Geometry
{
    // Geometry definition
    constexpr int NX = 64;
    constexpr int NY = 64;
    constexpr int NZ = 0;

    // Layer definition
    constexpr int LNX = NX;
    constexpr int LNY = 3;
    constexpr int LNZ = 0;
};

// Defining block size and number for kernel initialization, in this branch we use this part mainly for the initialization part.
#define BX 32
#define BY 16
#define GX (Geometry::NX / BX)
#define GY (Geometry::NY / BY)

inline dim3 block(BX, BY);
inline dim3 blockNumber(GX, GY);

// Defining block size and number for kernel using layer.
#define LBX Geometry::LNX
#define LBY 1
#define LGX (Geometry::LNX / LBX)
#define LGY (Geometry::LNY / LBY)

inline dim3 blockLayer = dim3(LBX, LBY);
inline dim3 blockNumberLayer = dim3(LGX, LGY);

// Stencil definition
namespace D2Q9
{
    constexpr int D = 2;
    constexpr int Q = 9;

    constexpr int momSize = Geometry::NX * Geometry::NY;
    constexpr int momByteSize = Geometry::NX * Geometry::NY * sizeof(float);

    constexpr int layerSize = Geometry::NX * Geometry::NY;
    constexpr int layerByteSize = Geometry::NX * Geometry::NY * sizeof(float);

    constexpr int momNum = 6; // Moment number for D2Q9 is set at 6. Look at momConfig.cuh for more details.

    const float a_s = sqrtf(3);
    constexpr float a_s2 = 3.f;
    constexpr float a_s4 = 9.f;
    constexpr float inv_as2 = 1.f / a_s2;

    __host__ __device__ __forceinline__ float w(int i)
    {
        switch (i)
        {
        case 0:
            return 4.f / 9.f;
        case 1:
        case 2:
        case 3:
        case 4:
            return 1.f / 9.f;
        case 5:
        case 6:
        case 7:
        case 8:
            return 1.f / 36.f;
        default:
            return 0.f;
        }
    }

    __host__ __device__ __forceinline__ float cx(int i)
    {
        switch (i)
        {
        case 0:
            return 0.f;
        case 1:
            return 1.f;
        case 2:
            return 0.f;
        case 3:
            return -1.f;
        case 4:
            return 0.f;
        case 5:
            return 1.f;
        case 6:
            return -1.f;
        case 7:
            return -1.f;
        case 8:
            return 1.f;
        default:
            return 0.f;
        }
    }

    __host__ __device__ __forceinline__ float cy(int i)
    {
        switch (i)
        {
        case 0:
            return 0.f;
        case 1:
            return 0.f;
        case 2:
            return 1.f;
        case 3:
            return 0.f;
        case 4:
            return -1.f;
        case 5:
            return 1.f;
        case 6:
            return 1.f;
        case 7:
            return -1.f;
        case 8:
            return -1.f;
        default:
            return 0.f;
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
