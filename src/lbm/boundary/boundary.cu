#include "boundary.cuh"

#include "../../config/mom_config.cuh"
#include "../build/build_grid.cuh"

#include "../core/grid_id.cuh"
#include "../core/from_id.cuh"
#include "../core/to_u8.cuh"
#include "../core/def_bc.cuh"

#include "../f_i.cuh"

__device__ void center(int x, int y,
                       D2Q9 sim)
{
    int index = grid_id();
    float rho = 0.f;
    float ux = 0.f;
    float uy = 0.f;
    float mxx = 0.f;
    float mxy = 0.f;
    float myy = 0.f;

#pragma unroll
    for (int i = 0; i < Q; i++)
    {
        int index_from = from_id(x, y, i);

        float fi = f_i(index_from, i, sim.layer);

        rho += fi;

        ux += fi * (float)c_ix[i];
        uy += fi * (float)c_iy[i];

        mxx += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxy += fi * (c_ix[i] * c_iy[i]);
        myy += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    sim.mom[momIdx<MomentId::rho>(index)] = rho;

    sim.mom[momIdx<MomentId::ux>(index)] = ux / rho;
    sim.mom[momIdx<MomentId::uy>(index)] = uy / rho;

    sim.mom[momIdx<MomentId::mxx>(index)] = mxx / rho;
    sim.mom[momIdx<MomentId::myy>(index)] = myy / rho;
    sim.mom[momIdx<MomentId::mxy>(index)] = mxy / rho;
}

__device__ void boundary(uint8_t mask_in, uint8_t node_in, int x, int y, D2Q9 sim)
{

    int index = grid_id();
    float sum_fi = 0.f;
    float mxy_I = 0.f;

    for (int i = 0; i < Q; i++)
    {
        if (mask_in & (1u << i))
        {
            int index_from = from_id(x, y, i);

            float fi = f_i(index_from, i, sim.layer);
            sum_fi += fi;
            mxy_I += fi * c_ix[i] * c_iy[i];
        }
    }

    mxy_I /= sum_fi;

    switch (static_cast<Boundary>(node_in))
    {
    case Boundary::East:
        Constants<Boundary::East> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        break;
    case Boundary::North:
        Constants<Boundary::North> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        break;
    case Boundary::West:
        Constants<Boundary::West> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        break;
    case Boundary::South:
        Constants<Boundary::South> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        break;
    case Boundary::Northeast:
        Constants<Boundary::Northeast> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        break;
    case Boundary::Northwest:
        Constants<Boundary::Northwest> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        break;
    case Boundary::Southwest:
        Constants<Boundary::Southwest> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        break;
    case Boundary::Southeast:
        Constants<Boundary::Southeast> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        break;
    default:
        break;
    }
}
