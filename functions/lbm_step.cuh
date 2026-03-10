#pragma once

__global__ void lbm_step(float* rhoA, float* uxA, float* uyA, float* mxxA, float* mxyA,  float* myyA,
                         float* rhoB, float* uxB, float* uyB, float* mxxB, float* mxyB,  float* myyB);