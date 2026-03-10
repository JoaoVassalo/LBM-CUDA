#pragma once

__global__ float propagation(int x, int y, float* rhoA, float* uxA, float* uyA, float* mxxA, float* mxyA,  float* myyA,
                                           float* rhoB, float* uxB, float* uyB, float* mxxB, float* mxyB,  float* myyB);