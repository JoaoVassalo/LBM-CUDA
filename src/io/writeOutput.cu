#include "writeOutput.cuh"

__global__ void packOutput(const varType *mom, varType *vtk)
{
    const size_t point = static_cast<size_t>(blockIdx.x) * blockDim.x + threadIdx.x;
    if (point >= D2Q9::momSize)
        return;

    const size_t count = D2Q9::momSize;
    const varType rho = mom[momIdx<momId::rho>(point)];
    const varType ux = mom[momIdx<momId::ux>(point)];
    const varType uy = mom[momIdx<momId::uy>(point)];

    vtk[point] = rho;
    vtk[count + point] = ux;
    vtk[2 * count + point] = uy;
    vtk[3 * count + 3 * point] = ux;
    vtk[3 * count + 3 * point + 1] = uy;
    vtk[3 * count + 3 * point + 2] = static_cast<varType>(0);
}

__host__ void writeOutput(moments sim, Grid2D grid, size_t t, size_t momByteSize, const std::string path)
{
    constexpr int blockSize = 256;
    const int blockCount = static_cast<int>((D2Q9::momSize + blockSize - 1) / blockSize);

    packOutput<<<blockCount, blockSize>>>(sim.mom, sim.vtk);
    cudaMemcpy(sim.mom_host, sim.vtk, momByteSize, cudaMemcpyDeviceToHost);
    writeVTI(t, path, sim.mom_host);

    if (t == 0)
    {
        return;
    }
}
