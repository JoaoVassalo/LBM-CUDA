#include "calc_density.cuh"
#include "calc_velocity.cuh"
#include "calc_momentum.cuh"
#include "collision.cuh"
#include "propagation.cuh"
#include "grid_id.cuh"


__global__ void lbm_step(float* rho, float* ux, float* uy, float* mxx, float* mxy,  float* myy){

    int x = blockIdx.x * blockDim.x + threadIdx.x;
    int y = blockIdx.y * blockDim.y + threadIdx.y;

    int index = grid_id(x, y);

    float mxx_col, mxy_col, myy_col = collision(index, rho, ux, uy, mxx, mxy, myy);

    



}