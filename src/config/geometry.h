#pragma once

struct Geometry
{ // Grid size.
    static const int Nx = 128;
    static const int Ny = 128;
    static const int grid_num = Nx * Ny;
};

// Simulation time.
constexpr int tf = 1e3;
constexpr int t_interval = tf / 100;