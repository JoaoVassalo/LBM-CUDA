#include "propagation.cuh"

#include <iostream>

#include "../presets/geometry.h"
#include "../presets/stencil.cuh"
#include "../presets/physics.h"

#include "pop_id.cuh"
#include "grid_id.cuh"
#include "from_id.cuh"
#include "./boundary/corner.cuh"
#include "./boundary/north.cuh"
#include "./boundary/wall.cuh"
#include "./boundary/center.cuh"
#include "./boundary/inlet.cuh"
#include "./boundary/outlet.cuh"
#include "./boundary/boundary.cuh"

__device__ void propagation(float *mom_in,
                            float *mom_out)
{
   int x = blockIdx.x * blockDim.x + threadIdx.x;
   int y = blockIdx.y * blockDim.y + threadIdx.y;

   int index = grid_id();

   if (x == 0 && y == 0)
   { // Sudoeste
      int size = 4;
      inlet_southwest(Is_SW, x, y,
                      mom_in, mom_out);
   }
   else if (x == Nx - 1 && y == 0)
   { // Sudeste
      int size = 4;
      wall_southeast(size, Is_SE, Os_SE, x, y,
                     mom_in,
                     mom_out);
   }
   else if (x == 0 && y == Ny - 1)
   { // Noroeste
      int size = 4;
      inlet_northwest(Is_NW, x, y,
                      mom_in,
                      mom_out);
   }
   else if (x == Nx - 1 && y == Ny - 1)
   { // Nordeste
      int size = 4;
      wall_northeast(size, Is_NE, Os_NE, x, y,
                     mom_in,
                     mom_out);
   }
   else if (y == 0)
   { // Sul
      int size = 6;
      wall_south(size, Is_S, Os_S, x, y,
                 mom_in,
                 mom_out);
   }
   else if (y == Ny - 1)
   { // Norte
      int size = 6;
      wall_north(size, Is_N, Os_N, x, y,
                 mom_in,
                 mom_out);
   }
   else if (x == 0)
   { // Oeste
      inlet_west(Is_W, x, y,
                 mom_in,
                 mom_out);
   }
   else if (x == Nx - 1)
   { // Leste
      int size = 6;
      wall_east(size, Is_E, Os_E, x, y,
                mom_in,
                mom_out);
   }
   else
   { // Centro

      center(x, y,
             mom_in,
             mom_out);
   }
}