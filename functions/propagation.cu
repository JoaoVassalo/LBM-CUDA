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

__device__ void propagation(varUnit *mom_in,
                            varUnit *mom_out)
{
   int x = blockIdx.x * blockDim.x + threadIdx.x;
   int y = blockIdx.y * blockDim.y + threadIdx.y;

   int index = grid_id();

   if (x == 0 && y == 0)
   { // Sudoeste
      int size = 4;
      // inlet_southwest(Is_SW, x, y, mom_in, mom_out);
      wall_southwest(size, Is_SW, Os_SW, x, y,
                     mom_in,
                     mom_out);
   }
   else if (x == Nx - 1 && y == 0)
   { // Sudeste
      int size = 4;
      // outlet_southeast(Is_SE, x, y, mom_in, mom_out);
      wall_southeast(size, Is_SE, Os_SE, x, y,
                     mom_in,
                     mom_out);
   }
   else if (x == 0 && y == Ny - 1)
   { // Noroeste
      int size = 4;
      // inlet_northwest(Is_SW, x, y, mom_in, mom_out);
      wall_northwest(size, Is_NW, Os_NW, x, y,
                     mom_in,
                     mom_out);
   }
   else if (x == Nx - 1 && y == Ny - 1)
   { // Nordeste
      int size = 4;
      // outlet_northeast(Is_NE, x, y, mom_in, mom_out);
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
      int size = 6;
      // wall_west(size, Is_W, Os_W, x, y, mom_in, mom_out);
      inlet_west(Is_W, x, y,
                 mom_in,
                 mom_out);
   }
   else if (x == Nx - 1)
   { // Leste
      int size = 6;
      outlet_east(Is_E, x, y,
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