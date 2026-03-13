#include "propagation.cuh"

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

__device__ void propagation(float *rho_in, float *ux_in, float *uy_in, float *mxx_in, float *mxy_in, float *myy_in,
                            float *rho_out, float *ux_out, float *uy_out, float *mxx_out, float *mxy_out, float *myy_out)
{
   int x = blockIdx.x * blockDim.x + threadIdx.x;
   int y = blockIdx.y * blockDim.y + threadIdx.y;

   int index = grid_id(x, y);

   if (x == 0 && y == 0)
   { // Sudoeste
      constexpr int I_s[4] = {0, 3, 4, 7};
      constexpr int O_s[4] = {0, 1, 2, 5};

      corner(I_s, O_s, x, y,
             rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in,
             rho_out);

      ux_out[index] = 0.f;
      uy_out[index] = 0.f;

      mxx_out[index] = 0.f;
      myy_out[index] = 0.f;
      mxy_out[index] = 0.f;
   }
   else if (x == Nx - 1 && y == 0)
   { // Sudeste
      constexpr int I_s[4] = {0, 1, 4, 8};
      constexpr int O_s[4] = {0, 2, 3, 6};

      corner(I_s, O_s, x, y,
             rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in,
             rho_out);

      ux_out[index] = 0.f;
      uy_out[index] = 0.f;

      mxx_out[index] = 0.f;
      myy_out[index] = 0.f;
      mxy_out[index] = 0.f;
   }
   else if (x == 0 && y == Ny - 1)
   { // Noroeste
      constexpr int size = 4;
      constexpr int I_s[size] = {0, 2, 3, 6};
      constexpr int O_s[size] = {0, 1, 4, 8};

      north(size, I_s, O_s, x, y,
            rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in,
            rho_out, mxy_out);

      ux_out[index] = u_max;
      uy_out[index] = 0;

      mxx_out[index] = u_max * u_max;
      myy_out[index] = 0.f;
   }
   else if (x == Nx - 1 && y == Ny - 1)
   { // Nordeste
      constexpr int size = 4;
      constexpr int I_s[size] = {0, 1, 2, 5};
      constexpr int O_s[size] = {0, 3, 4, 7};

      north(size, I_s, O_s, x, y,
            rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in,
            rho_out, mxy_out);

      ux_out[index] = u_max;
      uy_out[index] = 0.f;

      mxx_out[index] = u_max * u_max;
      myy_out[index] = 0.f;
   }
   else if (y == 0)
   { // Sul
      constexpr int I_s[6] = {0, 1, 3, 4, 7, 8};
      constexpr int O_s[6] = {0, 1, 2, 3, 5, 6};

      wall(I_s, O_s, x, y,
           rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in,
           rho_out, mxy_out);

      ux_out[index] = 0.f;
      uy_out[index] = 0.f;

      mxx_out[index] = 0.f;
      myy_out[index] = 0.f;
   }
   else if (y == Ny - 1)
   { // Norte
      constexpr int size = 6;
      constexpr int I_s[size] = {0, 1, 2, 3, 5, 6};
      constexpr int O_s[size] = {0, 1, 3, 4, 7, 8};

      north(size, I_s, O_s, x, y,
            rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in,
            rho_out, mxy_out);

      ux_out[index] = u_max;
      uy_out[index] = 0.f;

      mxx_out[index] = u_max * u_max;
      myy_out[index] = 0.f;
   }
   else if (x == 0)
   { // Oeste
      constexpr int I_s[6] = {0, 2, 3, 4, 6, 7};
      constexpr int O_s[6] = {0, 1, 2, 4, 5, 8};

      wall(I_s, O_s, x, y,
           rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in,
           rho_out, mxy_out);

      ux_out[index] = 0.f;
      uy_out[index] = 0.f;

      mxx_out[index] = 0.f;
      myy_out[index] = 0.f;
   }
   else if (x == Nx - 1)
   { // Leste
      constexpr int I_s[6] = {0, 1, 2, 4, 5, 8};
      constexpr int O_s[6] = {0, 2, 3, 4, 6, 7};

      wall(I_s, O_s, x, y,
           rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in,
           rho_out, mxy_out);

      ux_out[index] = 0.f;
      uy_out[index] = 0.f;

      mxx_out[index] = 0.f;
      myy_out[index] = 0.f;
   }
   else
   { // Centro
      center(x, y,
             rho_in, ux_in, uy_in, mxx_in, mxy_in, myy_in,
             rho_out, ux_out, uy_out, mxx_out, mxy_out, myy_out);
   }
}