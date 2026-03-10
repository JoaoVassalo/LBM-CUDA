#include "calc_density.cuh"
#include "calc_velocity.cuh"
#include "calc_momentum.cuh"
#include "collision.cuh"
#include "propagation.cuh"
#include "grid_id.cuh"


__global__ void lbm_step(float* rho_in, float* ux_in, float* uy_in, float* mxx_in, float* mxy_in,  float* myy_in,
                         float* rho_out, float* ux_out, float* uy_out, float* mxx_out, float* mxy_out,  float* myy_out){

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    int index = grid_id(x, y);

    collision(index, ux_in, uy_in, mxx_in, mxy_in, myy_in);

    propagation(x, y, rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in, 
                      rho_out, ux_out, uy_out, mxx_out, mxy_out, myy_out);

    
    
}