#pragma once

// Grid size.
constexpr int Nx = 128;
constexpr int Ny = 128;
constexpr int grid_num = Nx * Ny;

// Simulation time.
constexpr int tf = 1e3;
constexpr int t_interval = tf / 100;

// Block size
#define BX 32
#define BY 16
#define GX (Nx / BX)
#define GY (Ny / BY)