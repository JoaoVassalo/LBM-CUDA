#include "writeOutput.cuh"

__host__ void writeOutput(moments sim, Grid2D grid, size_t t, size_t momByteSize, const std::string path)
{
    cudaDeviceSynchronize();
    cudaMemcpy(sim.mom_host, sim.mom, momByteSize, cudaMemcpyDeviceToHost);
    writeVTI(t, path, sim.mom_host);

    std::cout << "Iteration " << t << std::endl;
    std::cout << std::endl;

    if (t == 0)
    {
        return;
    }
}