#pragma once

#include "presets/stencil.cuh"
#include "presets/physics.h"

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

template <Boundary bc>
struct Constants
{
public:
    __host__ void calc_constant(float ux, float uy, int size)
    {
#pragma unroll
        for (int k = 0; k < size; k++)
        {
            // Os
            int i = Os_E[k];

            C1 += w[i] *
                      (1.f +
                       a_s2 * ux * c_ix[i] +
                       a_s2 * uy * c_iy[i] +
                       a_s4 * 0.5f * ux * ux +
                       a_s4 * 0.5f * uy * uy) +
                  omega * a_s4 * 0.5f * ux * uy * c_ix[i] * c_iy[i];
            C2 += (1.f - omega) * w[i] * a_s4 * 0.5f * c_ix[i] * c_iy[i];

            /*--------------------------------------------------------------------------*/
            // Is
            i = Is_E[k];

            C3 += w[i] *
                  (1.f +
                   a_s2 * ux * c_ix[i] + a_s2 * uy * c_iy[i] +
                   a_s4 * 0.5f * ux * ux * (c_ix[i] * c_ix[i] - inv_as2) +
                   a_s4 * 0.5f * uy * uy * (c_iy[i] * c_iy[i] - inv_as2)) *
                  c_ix[i] * c_iy[i];
            C4 += w[i] * a_s4 * c_ix[i] * c_iy[i] * c_ix[i] * c_iy[i];
        }

        return;
    }
};

// East
template <>
struct Constants<Boundary::East>
{
    static constexpr int size = 6;

    static constexpr float ux = 0.f;
    static constexpr float uy = 0.f;

    calc_constant(ux, uy);
};
// North
template <>
struct Constants<Boundary::North>
{
    static constexpr int size = 6;

    static constexpr float ux = 1.f;
    static constexpr float uy = 0.f;

    calc_constant(ux, uy);
};
// West
template <>
struct Constants<Boundary::West>
{
    static constexpr int size = 6;

    static constexpr float ux = 0.f;
    static constexpr float uy = 0.f;

    calc_constant(ux, uy);
};
// South
template <>
struct Constants<Boundary::South>
{
    static constexpr int size = 6;

    static constexpr float ux = 0.f;
    static constexpr float uy = 0.f;

    calc_constant(ux, uy);
};
// Northeast
template <>
struct Constants<Boundary::Northeast>
{
    static constexpr int size = 4;

    static constexpr float ux = 1.f;
    static constexpr float uy = 0.f;

    calc_constant(ux, uy);
};
// Northwest
template <>
struct Constants<Boundary::Northwest>
{
    static constexpr int size = 4;

    static constexpr float ux = 1.f;
    static constexpr float uy = 0.f;

    calc_constant(ux, uy);
};
// Southwest
template <>
struct Constants<Boundary::Southwest>
{
    static constexpr int size = 4;

    static constexpr float ux = 0.f;
    static constexpr float uy = 0.f;

    calc_constant(ux, uy);
};
// Southeast
template <>
struct Constants<Boundary::Southeast>
{
    static constexpr int size = 4;

    static constexpr float ux = 0.f;
    static constexpr float uy = 0.f;

    calc_constant(ux, uy);
};
