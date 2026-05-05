#pragma once

// Grid size.
constexpr int Nx = 64;
constexpr int Ny = 64;
constexpr int grid_num = Nx * Ny;

// Layer size.
constexpr int LNx = Nx;
constexpr int LNy = 3;
constexpr int layer_num = LNx * LNy;

// Simulation time.
constexpr int tf = 4e6;
constexpr int t_interval = 4e4;

// Block size
#define BX 32
#define BY 16
#define GX (Nx / BX)
#define GY (Ny / BY)