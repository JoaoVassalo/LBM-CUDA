#pragma once

#include <math.h>

struct Geometry
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

// Stencil definition
struct D2Q9
{
    static const int D = 2;
    static const int Q = 9;

    float *mom;

    float *layer[Geometry::LNY];

    static const int momSize = Geometry::NX * Geometry::NY;
    static const int momByteSize = Geometry::NX * Geometry::NY * sizeof(float);

    static const int layerSize = Geometry::NX * Geometry::NY;
    static const int layerByteSize = Geometry::NX * Geometry::NY * sizeof(float);

    static const int momNum = 6; // Moment number for D2Q9 is set at 6. Look at momConfig.cuh for more details.

    __constant__ float w[Q];
    __constant__ float c_ix[Q];
    __constant__ float c_iy[Q];

    const float a_s = sqrtf(3);
    const float a_s2 = 3.f;
    const float a_s4 = 9.f;
    const float inv_as2 = 1.f / a_s2;
};
