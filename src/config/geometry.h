#pragma once

struct Geometry
{ // Grid size.
    static const int Nx = 64;
    static const int Ny = 64;
    static const int grid_num = Nx * Ny;
};

// Simulation time.
constexpr int tf = 4e6;
constexpr int t_interval = 4e4;
