#include "propagation.cuh"

__global__ void propagation(D2Q9 sim, layer layer, Grid2D grid)
{
   const int x = blockIdx.x * blockDim.x + threadIdx.x;
   const int y = blockIdx.y * blockDim.y + threadIdx.y;

   const int index = grid_id();

   if (grid.node[index] & to_u8(domainTags::Boundary))
   {
      boundary(grid.mask[index], grid.node[index], x, y, sim, layer);
   }
   else
   {
      center(x, y, sim, layer);
   }
}