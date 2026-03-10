#include "../presets/geometry.h"
#include "../presets/physics.h"

#include "collision.cuh"


__device__ float collision(int index, float* rho, float* ux, float* uy, float* mxx, float* mxy,  float* myy){
    
    if (index>Nx*Ny) return;

    float mxx_col = ((1-omega)*mxx[index] - omega*ux[index]*ux[index]);
    float mxy_col = ((1-omega)*mxy[index] - omega*ux[index]*uy[index]);
    float myy_col = ((1-omega)*myy[index] - omega*uy[index]*uy[index]);

    return mxx_col, mxy_col, myy_col;

}