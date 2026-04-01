__host__ void memory_init(int size)
{
    float *rhoA,
        *rhoB;
    float *uxA;
    float *uxB;
    float *uyA;
    float *uyB;
    float *mxxA;
    float *mxxB;
    float *mxyA;
    float *mxyB;
    float *myyA;
    float *myyB;

    cudaMalloc((void **)&rhoA, size);
    cudaMalloc((void **)&rhoB, size);
    cudaMalloc((void **)&uxA, size);
    cudaMalloc((void **)&uxB, size);
    cudaMalloc((void **)&uyA, size);
    cudaMalloc((void **)&uyB, size);
    cudaMalloc((void **)&mxxA, size);
    cudaMalloc((void **)&mxxB, size);
    cudaMalloc((void **)&mxyA, size);
    cudaMalloc((void **)&mxyB, size);
    cudaMalloc((void **)&myyA, size);
    cudaMalloc((void **)&myyB, size);
}