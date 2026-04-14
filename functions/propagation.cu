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

__device__ void propagation(float *mom_in,
                            float *mom_out)
{
   int x = blockIdx.x * blockDim.x + threadIdx.x;
   int y = blockIdx.y * blockDim.y + threadIdx.y;

   int index = grid_id();

   if (x == 0 && y == 0)
   { // Sudoeste

      corner(Is_SW, Os_SW, x, y,
             mom_in,
             mom_out);

      mom_out[momIdx<MomentId::ux>(index)] = 0.f; // ux
      mom_out[momIdx<MomentId::uy>(index)] = 0.f; // uy

      mom_out[momIdx<MomentId::mxx>(index)] = 0.f; // mxx
      mom_out[momIdx<MomentId::myy>(index)] = 0.f; // myy
      mom_out[momIdx<MomentId::myy>(index)] = 0.f; // mxy
   }
   else if (x == Nx - 1 && y == 0)
   { // Sudeste

      corner(Is_SE, Os_SE, x, y,
             mom_in,
             mom_out);

      mom_out[momIdx<MomentId::ux>(index)] = 0.f; // ux
      mom_out[momIdx<MomentId::uy>(index)] = 0.f; // uy

      mom_out[momIdx<MomentId::mxx>(index)] = 0.f; // mxx
      mom_out[momIdx<MomentId::myy>(index)] = 0.f; // myy
      mom_out[momIdx<MomentId::mxy>(index)] = 0.f; // mxy
   }
   else if (x == 0 && y == Ny - 1)
   { // Noroeste
      constexpr int size = 4;

      north(size, Is_NW, Os_NW, x, y,
            mom_in,
            mom_out);

      mom_out[momIdx<MomentId::ux>(index)] = u_max; // ux
      mom_out[momIdx<MomentId::uy>(index)] = 0;     // uy

      mom_out[momIdx<MomentId::mxx>(index)] = u_max * u_max; // mxx
      mom_out[momIdx<MomentId::myy>(index)] = 0.f;           // myy
   }
   else if (x == Nx - 1 && y == Ny - 1)
   { // Nordeste
      constexpr int size = 4;

      north(size, Is_NE, Os_NE, x, y,
            mom_in,
            mom_out);

      mom_out[momIdx<MomentId::ux>(index)] = u_max; // ux
      mom_out[momIdx<MomentId::uy>(index)] = 0.f;   // uy

      mom_out[momIdx<MomentId::mxx>(index)] = u_max * u_max; // mxx
      mom_out[momIdx<MomentId::myy>(index)] = 0.f;           // myy
   }
   else if (y == 0)
   { // Sul

      wall(Is_S, Os_S, x, y,
           mom_in,
           mom_out);

      mom_out[momIdx<MomentId::ux>(index)] = 0.f; // ux
      mom_out[momIdx<MomentId::uy>(index)] = 0.f; // uy

      mom_out[momIdx<MomentId::mxx>(index)] = 0.f; // mxx
      mom_out[momIdx<MomentId::myy>(index)] = 0.f; // myy
   }
   else if (y == Ny - 1)
   { // Norte
      constexpr int size = 6;

      north(size, Is_N, Os_N, x, y,
            mom_in,
            mom_out);

      mom_out[momIdx<MomentId::ux>(index)] = u_max; // ux
      mom_out[momIdx<MomentId::uy>(index)] = 0.f;   // uy

      mom_out[momIdx<MomentId::mxx>(index)] = u_max * u_max; // mxx
      mom_out[momIdx<MomentId::myy>(index)] = 0.f;           // myy
   }
   else if (x == 0)
   { // Oeste

      wall(Is_W, Os_W, x, y,
           mom_in,
           mom_out);

      mom_out[momIdx<MomentId::ux>(index)] = 0.f; // ux
      mom_out[momIdx<MomentId::uy>(index)] = 0.f; // uy

      mom_out[momIdx<MomentId::mxx>(index)] = 0.f; // mxx
      mom_out[momIdx<MomentId::myy>(index)] = 0.f; // myy
   }
   else if (x == Nx - 1)
   { // Leste

      wall(Is_E, Os_E, x, y,
           mom_in,
           mom_out);

      mom_out[momIdx<MomentId::ux>(index)] = 0.f; // ux
      mom_out[momIdx<MomentId::uy>(index)] = 0.f; // uy

      mom_out[momIdx<MomentId::mxx>(index)] = 0.f; // mxx
      mom_out[momIdx<MomentId::myy>(index)] = 0.f; // myy
   }
   else
   { // Centro

      center(x, y,
             mom_in,
             mom_out);
   }
}