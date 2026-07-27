#pragma once

// Grid size.
constexpr int Nx = 512;
constexpr int Ny = 128;
constexpr int grid_num = Nx * Ny;

// Simulation time.
constexpr int tf = 1e2;
constexpr int printNum = 100;
constexpr int t_interval = tf / printNum;