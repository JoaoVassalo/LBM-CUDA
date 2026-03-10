#include <iostream>

#include "presets/stencil.h"
#include "presets/geometry.h"
#include "presets/config.h"

#include "functions/init_domain.cuh"
#include "functions/lbm_step.cuh"


int main(){

    dim3 block(Nx, Ny);

    float *rhoA; float *rhoB;
    float *uxA; float *uxB; 
    float *uyA; float *uyB; 
    float *mxxA; float *mxxB; 
    float *mxyA; float *mxyB; 
    float *myyA; float *myyB; 

    int size = Nx*Ny*sizeof(float);

    

    cudaMalloc((void**)&rhoA, size); cudaMalloc((void**)&rhoB, size);
    cudaMalloc((void**)&uxA, size); cudaMalloc((void**)&uxB, size);
    cudaMalloc((void**)&uyA, size); cudaMalloc((void**)&uyB, size);
    cudaMalloc((void**)&mxxA, size); cudaMalloc((void**)&mxxB, size);
    cudaMalloc((void**)&mxyA, size); cudaMalloc((void**)&mxyB, size);
    cudaMalloc((void**)&myyA, size); cudaMalloc((void**)&myyB, size);
    
    initDomain<<<N_block, block>>>(rhoA, uxA, uyA, mxxA, mxyA, myyA);

    cudaError_t err = cudaGetLastError();
    printf("Kernel launch error: %s\n", cudaGetErrorString(err));

    cudaDeviceSynchronize();

    cudaMemcpy(rhoB, rhoA, size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(uxB, uxA, size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(uyB, uyA, size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(mxxB, mxxA, size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(mxyB, mxyA, size, cudaMemcpyDeviceToDevice);
    cudaMemcpy(myyB, myyA, size, cudaMemcpyDeviceToDevice);

    for (int t = 0; t<tf; t++){
        bool id = (t & 1);
        
        if (!id){
            lbm_step<<<N_block, block>>>(rhoA, uxA, uyA, mxxA, mxyA, myyA);
        }
        else{
            lbm_step<<<N_block, block>>>(rhoB, uxB, uyB, mxxB, mxyB, myyB);
        }

        cudaDeviceSynchronize();
    }

}