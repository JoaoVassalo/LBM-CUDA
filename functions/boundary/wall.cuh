#pragma once

#include "../../presets/config.h"

__device__ void wall(CInt *I_s, CInt *O_s, int x, int y,
                     varUnit *mom_in,
                     varUnit *mom_out);

__device__ void wall_north(int size, CInt *I_s, CInt *O_s,
                           int x, int y,
                           varUnit *mom_in,
                           varUnit *mom_out);

__device__ void wall_south(int size, CInt *I_s, CInt *O_s,
                           int x, int y,
                           varUnit *mom_in,
                           varUnit *mom_out);

__device__ void wall_east(int size, CInt *I_s, CInt *O_s,
                          int x, int y,
                          varUnit *mom_in,
                          varUnit *mom_out);

__device__ void wall_west(int size, CInt *I_s, CInt *O_s,
                          int x, int y,
                          varUnit *mom_in,
                          varUnit *mom_out);

__device__ void wall_northeast(int size, CInt *I_s, CInt *O_s,
                               int x, int y,
                               varUnit *mom_in,
                               varUnit *mom_out);

__device__ void wall_northwest(int size, CInt *I_s, CInt *O_s,
                               int x, int y,
                               varUnit *mom_in,
                               varUnit *mom_out);

__device__ void wall_southeast(int size, CInt *I_s, CInt *O_s,
                               int x, int y,
                               varUnit *mom_in,
                               varUnit *mom_out);

__device__ void wall_southwest(int size, CInt *I_s, CInt *O_s,
                               int x, int y,
                               varUnit *mom_in,
                               varUnit *mom_out);