#include <propagation.cuh>

#include "../presets/geometry.h"
#include "../presets/stencil.h"


#include <pop_id.cuh>
#include <grid_id.cuh>
#include <from_id.cuh>
#include <calc_density.cuh>


__device__ void propagation(float* rho_in, float* ux_in, float* uy_in, float* mxx_in, float* mxy_in,  float* myy_in,
                             float* rho_out, float* ux_out, float* uy_out, float* mxx_out, float* mxy_out,  float* myy_out){
 


   int x = blockIdx.x * blockDim.x + threadIdx.x;
   int y = blockIdx.y * blockDim.y + threadIdx.y;
   
   if (x >= Nx || y >= Ny) return;

   int index = grid_id(x,y);

   rho_out[index] = 0;
   ux_out[index] = 0;
   uy_out[index] = 0;
   mxx_out[index] = 0;
   mxy_out[index] = 0;
   myy_out[index] = 0;


   #pragma unroll
   for (int i = 0; i < Q; i++) {
      int index_from = from_id(x, y, i);
      
      float f = rho_in[index_from]*w[i]*(1 + a_s2*ux_in[index_from]*c_ix[i] + a_s2*uy_in[index_from]*c_iy[i] + 
                                                                  a_s4*0.5f*mxx_in[index_from]*(c_ix[i]*c_ix[i] - inv_as2) + 
                                                                  a_s4*0.5f*mxy_in[index_from]*(c_ix[i]*c_iy[i]) + 
                                                                  a_s4*0.5f*myy_in[index_from]*(c_iy[i]*c_iy[i] - inv_as2) );

      rho_out[index] += f;

      ux_out[index] += c_ix[i]*f;

      uy_out[index] += c_iy[i]*f;

      mxx_out[index] += (c_ix[i]*c_ix[i] - inv_as2)*f;

      mxy_out[index] += (c_ix[i]*c_iy[i])*f;

      myy_out[index] += (c_iy[i]*c_iy[i] - inv_as2)*f;


   }

   

}