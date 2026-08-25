# LBM-CUDA

My first Lattice Boltzmann Method (LBM) implementation using **CUDA and C++**. This project was developed as a first step toward GPU-accelerated LBM simulations, using a moment-based formulation and the D2Q9 lattice.

## Method

The solver uses the **D2Q9** lattice and stores the macroscopic and non-equilibrium information in a moment-based representation:

rho, ux, uy, mxx, mxy, myy

The main simulation loop performs propagation and collision on the GPU. The collision step uses a BGK relaxation parameter, while the propagation step handles the lattice transport and the cavity boundary conditions.

The current configuration is set up for a **2D square cavity with a moving lid**.

## CUDA implementation

The simulation is parallelized using CUDA kernels. The domain is divided into CUDA blocks, with the default configuration using:

- `BX = 32`
- `BY = 16`
- `Nx = 1024`
- `Ny = 1024`

The solver uses two GPU memory buffers (`momA` and `momB`) for the time evolution.

Performance is measured in **MLUPS (Million Lattice Updates Per Second)**.

## Structure

- `main.cu` — main simulation loop, memory allocation, execution, output, and performance measurement.
- `functions/` — LBM operations, including initialization, propagation, collision, indexing, and boundary conditions.
- `functions/boundary/` — cavity wall and moving-lid boundary implementations.
- `functions/streaming/` — streaming-related routines.
- `presets/` — geometry, physical parameters, D2Q9 stencil, and CUDA configuration.
- `cuda_config/` — simulation data structures and configuration.
- `animation/` — Python scripts for plotting and generating animations.
- `vtk.cu/.cuh` — VTK output for visualization.
- `run.sh` — compilation and execution script.

## Execution

CUDA Toolkit and `nvcc` are required.

Run the simulation with:

```bash
bash run.sh