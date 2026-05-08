#include "boundary.cuh"

#include "../../config/mom_config.cuh"
#include "../build/build_grid.cuh"

#include "../core/grid_id.cuh"
#include "../core/to_u8.cuh"

__device__ void applyBoundary(uint8_t *node, uint8_t *mask)
{
    int index = grid_id();
}
