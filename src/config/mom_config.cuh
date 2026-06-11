#pragma once

#include <cstdint>
#include <stdio.h>

#include "stencil.cuh"
#include "physics.h"
#include "../lbm/build/build_grid.cuh"

enum Boundary
{
    Center = 0,
    East,
    North,
    West,
    South,
    Northeast,
    Northwest,
    Southwest,
    Southeast
};

__device__ __forceinline__ void calc_constant(float ux, float uy, uint8_t mask,
                                              float &C1, float &C2, float &C3, float &C4)
{
    C1 += w[0] *
              (1.f +
               a_s2 * ux * c_ix[0] +
               a_s2 * uy * c_iy[0] +
               a_s4 * 0.5f * ux * ux +
               a_s4 * 0.5f * uy * uy) +
          omega * a_s4 * 0.5f * ux * uy * c_ix[0] * c_iy[0];
    C2 += (1.f - omega) * w[0] * a_s4 * 0.5f * c_ix[0] * c_iy[0];
    C3 += w[0] *
          (1.f +
           a_s2 * ux * c_ix[0] + a_s2 * uy * c_iy[0] +
           a_s4 * 0.5f * ux * ux * (c_ix[0] * c_ix[0] - inv_as2) +
           a_s4 * 0.5f * uy * uy * (c_iy[0] * c_iy[0] - inv_as2)) *
          c_ix[0] * c_iy[0];
    C4 += w[0] * a_s4 * c_ix[0] * c_iy[0] * c_ix[0] * c_iy[0];
    for (int k = 1; k < Q; k++)
    {
        if (mask & (1u << (k - 1)))
        {
            // Os
            C1 += w[k] *
                      (1.f +
                       a_s2 * ux * c_ix[k] +
                       a_s2 * uy * c_iy[k] +
                       a_s4 * 0.5f * ux * ux +
                       a_s4 * 0.5f * uy * uy) +
                  omega * a_s4 * 0.5f * ux * uy * c_ix[k] * c_iy[k];
            C2 += (1.f - omega) * w[k] * a_s4 * 0.5f * c_ix[k] * c_iy[k];

            /*--------------------------------------------------------------------------*/
            // Is
            int i = income(k);

            C3 += w[i] *
                  (1.f +
                   a_s2 * ux * c_ix[i] + a_s2 * uy * c_iy[i] +
                   a_s4 * 0.5f * ux * ux * (c_ix[i] * c_ix[i] - inv_as2) +
                   a_s4 * 0.5f * uy * uy * (c_iy[i] * c_iy[i] - inv_as2)) *
                  c_ix[i] * c_iy[i];
            C4 += w[i] * a_s4 * c_ix[i] * c_iy[i] * c_ix[i] * c_iy[i];
        }
    }
}

template <Boundary bc>
struct Constants
{
};

// East
template <>
struct Constants<Boundary::East>
{

    static constexpr float ux = 0.f;
    static constexpr float uy = 0.f;

    float C1 = 0.f,
          C2 = 0.f,
          C3 = 0.f,
          C4 = 0.f;
    __device__ Constants(uint8_t mask)
    {
        calc_constant(ux, uy, mask, C1, C2, C3, C4);
    }
};

// North
template <>
struct Constants<Boundary::North>
{

    static constexpr float ux = u_max;
    static constexpr float uy = 0.f;

    float C1 = 0.f,
          C2 = 0.f,
          C3 = 0.f,
          C4 = 0.f;

    __device__ Constants(uint8_t mask)
    {
        calc_constant(ux, uy, mask, C1, C2, C3, C4);
    }
};
// West
template <>
struct Constants<Boundary::West>
{

    static constexpr float ux = 0.f;
    static constexpr float uy = 0.f;

    float C1 = 0.f,
          C2 = 0.f,
          C3 = 0.f,
          C4 = 0.f;

    __device__ Constants(uint8_t mask)
    {
        calc_constant(ux, uy, mask, C1, C2, C3, C4);
    }
};
// South
template <>
struct Constants<Boundary::South>
{

    static constexpr float ux = 0.f;
    static constexpr float uy = 0.f;

    float C1 = 0.f,
          C2 = 0.f,
          C3 = 0.f,
          C4 = 0.f;

    __device__ Constants(uint8_t mask)
    {
        calc_constant(ux, uy, mask, C1, C2, C3, C4);
    }
};
// Northeast
template <>
struct Constants<Boundary::Northeast>
{

    static constexpr float ux = u_max;
    static constexpr float uy = 0.f;

    float C1 = 0.f,
          C2 = 0.f,
          C3 = 0.f,
          C4 = 0.f;

    __device__ Constants(uint8_t mask)
    {
        calc_constant(ux, uy, mask, C1, C2, C3, C4);
    }
};
// Northwest
template <>
struct Constants<Boundary::Northwest>
{

    static constexpr float ux = u_max;
    static constexpr float uy = 0.f;

    float C1 = 0.f,
          C2 = 0.f,
          C3 = 0.f,
          C4 = 0.f;

    __device__ Constants(uint8_t mask)
    {
        calc_constant(ux, uy, mask, C1, C2, C3, C4);
    }
};
// Southwest
template <>
struct Constants<Boundary::Southwest>
{

    static constexpr float ux = 0.f;
    static constexpr float uy = 0.f;

    float C1 = 0.f,
          C2 = 0.f,
          C3 = 0.f,
          C4 = 0.f;

    __device__ Constants(uint8_t mask)
    {
        calc_constant(ux, uy, mask, C1, C2, C3, C4);
    }
};
// Southeast
template <>
struct Constants<Boundary::Southeast>
{

    static constexpr float ux = 0.f;
    static constexpr float uy = 0.f;

    float C1 = 0.f,
          C2 = 0.f,
          C3 = 0.f,
          C4 = 0.f;

    __device__ Constants(uint8_t mask)
    {
        calc_constant(ux, uy, mask, C1, C2, C3, C4);
    }
};
