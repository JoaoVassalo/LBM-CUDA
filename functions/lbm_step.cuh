#pragma once

__global__ void lbm_step(float* rho, float* ux, float* uy, float* mxx, float* mxy,  float* myy);