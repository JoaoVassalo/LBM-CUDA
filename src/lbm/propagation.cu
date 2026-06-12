#include "propagation.cuh"

__global__ void propagate_layer_at(D2Q9 sim, layer current_layer, Grid2D grid, int y)
{
   const int x = blockIdx.x * blockDim.x + threadIdx.x;

   if (x >= Geometry::Nx)
      return;

   const int index = grid_id(x, y);

   if (grid.node[index] != to_u8(domainTags::Boundary))
   {
      boundary(grid.mask[index], grid.node[index], x, y, index, sim, current_layer);
   }
   else
   {
      center(grid.mask[index], x, y, index, sim, current_layer);
   }
}