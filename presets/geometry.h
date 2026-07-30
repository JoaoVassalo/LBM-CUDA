#pragma once

// Grid size.
constexpr int Nx = 512;
constexpr int Ny = 128;
constexpr int grid_num = Nx * Ny;

// Simulation time.
constexpr int tf = 4e6;
constexpr int printNum = 10;
constexpr int t_interval = tf / printNum;