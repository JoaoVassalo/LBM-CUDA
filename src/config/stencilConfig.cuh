#pragma once

#include <math.h>

namespace Geometry
{
    // Geometry definition
    static const int NX = 64;
    static const int NY = 64;
    static const int NZ = 0;

    // Layer definition
    static const int LNX = NX;
    static const int LNY = 3;
    static const int LNZ = 0;
};

// Defining block size and number for kernel initialization, in this branch we use this part mainly for the initialization part.
#define BX 32
#define BY 16
#define GX (Geometry::NX / BX)
#define GY (Geometry::NY / BY)

dim3 block = dim3(BX, BY);
dim3 blockNumber = dim3(GX, GY);

// Defining block size and number for kernel using layer.
#define LBX Geometry::LNX
#define LBY 1
#define LGX (Geometry::LNX / LBX)
#define LGY (Geometry::LNY / LBY)

dim3 blockLayer = dim3(LBX, LBY);
dim3 blockNumberLayer = dim3(LGX, LGY);

// Stencil definition
namespace D2Q9
{
    static const int D = 2;
    static const int Q = 9;

    static const int momSize = Geometry::NX * Geometry::NY;
    static const int momByteSize = Geometry::NX * Geometry::NY * sizeof(float);

    static const int layerSize = Geometry::NX * Geometry::NY;
    static const int layerByteSize = Geometry::NX * Geometry::NY * sizeof(float);

    static const int momNum = 6; // Moment number for D2Q9 is set at 6. Look at momConfig.cuh for more details.

    const float a_s = sqrtf(3);
    const float a_s2 = 3.f;
    const float a_s4 = 9.f;
    const float inv_as2 = 1.f / a_s2;

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
        }
    }
};
