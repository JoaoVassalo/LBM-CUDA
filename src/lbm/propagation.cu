#include "propagation.cuh"

#include <iostream>

#include "../config/geometry.h"
#include "../config/stencil.cuh"
#include "../config/physics.h"

#include "build/build_grid.cuh"

#include "../core/to_u8.cuh"

#include "boundary/boundary.cuh"

#include "pop_id.cuh"
#include "grid_id.cuh"
#include "from_id.cuh"

__global__ void propagation(D2Q9 sim, Grid2D grid)
{
   const int x = blockIdx.x * blockDim.x + threadIdx.x;
   const int y = blockIdx.y * blockDim.y + threadIdx.y;

   const int index = grid_id();

   if (grid.node[index] & to_u8(domainTags::Boundary))
   {
      boundary(grid.mask[index], grid.node[index], x, y, sim);
   }
   else
   {
      center(x, y, sim);
   }
}