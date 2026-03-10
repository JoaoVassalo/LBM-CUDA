#include <propagation.cuh>

#include <pop_id.cuh>
#include <grid_id.cuh>
#include "../presets/geometry.h"
#include "../presets/stencil.h"


__device__ float propagation(float x, float y, float rho, float ux, float uy, float mxx, float mxy,  float myy){
 


   int x = blockIdx.x + blockDim.x + threadIdx.x;
   int y = blockIdx.y + blockDim.y + threadIdx.y;
   
   if (x >= Nx || y >= Ny) return;


   #pragma unroll
   for (int i = 0; i < Q; i++) {
      
      int x_to = (x + c_ix[i] + Nx)%Nx;
      int y_to = (y + c_iy[i] + Ny)%Ny;

      int idx_to = pop_id(grid_id(x_to, y_to), i);
      int idx_from = pop_id(grid_id(x, y), i);

      f_i[idx_to] = f_col[idx_from];

   }

   

}