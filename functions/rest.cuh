#pragma once

__host__ __device__ __forceinline__ int rest(int a, int b)
{
    return a - b * (a / b);
}