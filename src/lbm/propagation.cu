#include "propagation.cuh"

#include <iostream>

#include "../config/geometry.h"
#include "../config/stencil.cuh"
#include "../config/physics.h"

#include "build/build_grid.cuh"

#include "../core/to_u8.cuh"

#include "boundary/boundary.cuh"

#include "pop_id.cuh"
#include "grid_id.cuh"
#include "from_id.cuh"

__device__ void propagation(float *mom, uint8_t *mask, uint8_t *node)
{
   int index = grid_id();

   if (node[index] & to_u8(domainTags::Fluid))
   {
      applyBoundary(mask);
   }
}