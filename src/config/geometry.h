#pragma once

struct Geometry
{ // Grid size.
    static const int Nx = 64;
    static const int Ny = 64;
    static const int grid_num = Nx * Ny;
};

// Block size
#define BX 32
#define BY 16
#define GX (Geometry::Nx / BX)
#define GY (Geometry::Ny / BY)

// Layer block size
#define LBX 256
#define LBY 1
#define LGX (Geometry::Nx / LBX)
#define LGY (Geometry::Ny / LBY)

struct layer
{
    static const int LNx = Geometry::Nx;
    static const int LNy = 3;
    static const int layer_num = Geometry::Ny / LNy;
};

// Simulation time.
constexpr int tf = 4e6;
constexpr int t_interval = 4e4;
