#include "boundary.cuh"

__device__ void center(uint8_t mask_in, int x, int y, int index,
                       D2Q9 sim, layer current_layer)
{
    float rho = 0.f;
    float ux = 0.f;
    float uy = 0.f;
    float mxx = 0.f;
    float mxy = 0.f;
    float myy = 0.f;

    float fi = f_i<0>(x, current_layer);

    rho += fi;

    ux += fi * (float)c_ix[0];
    uy += fi * (float)c_iy[0];

    mxx += fi * (c_ix[0] * c_ix[0] - inv_as2);
    mxy += fi * (c_ix[0] * c_iy[0]);
    myy += fi * (c_iy[0] * c_iy[0] - inv_as2);

    constexpr_for<1, Q>([&](auto qi)
                        {
        constexpr int i = decltype(qi)::value;

        if (mask_in & (1u << (i - 1)))
        {
            fi = f_i<i>(x, current_layer);

            rho += fi;

            ux += fi * (float)c_ix[i];
            uy += fi * (float)c_iy[i];

            mxx += fi * (c_ix[i] * c_ix[i] - inv_as2);
            mxy += fi * (c_ix[i] * c_iy[i]);
            myy += fi * (c_iy[i] * c_iy[i] - inv_as2);
        } });
    sim.mom[momIdx<MomentId::rho>(index)] = rho;

    sim.mom[momIdx<MomentId::ux>(index)] = ux / rho;
    sim.mom[momIdx<MomentId::uy>(index)] = uy / rho;

    sim.mom[momIdx<MomentId::mxx>(index)] = mxx / rho;
    sim.mom[momIdx<MomentId::myy>(index)] = myy / rho;
    sim.mom[momIdx<MomentId::mxy>(index)] = mxy / rho;
}

__device__ void boundary(uint8_t mask_in, uint8_t node_in,
                         int x, int y, int index,
                         D2Q9 sim, layer current_layer)
{
    float sum_fi = 0.f;
    float mxy_I = 0.f;

    float fi = f_i<0>(x, current_layer);
    sum_fi += fi;
    mxy_I += fi * c_ix[0] * c_iy[0];

    constexpr_for<1, Q>([&](auto qi)
                        {
        constexpr int i = decltype(qi)::value;

        if (mask_in & (1u << (i - 1)))
        {

            fi = f_i<i>(x, current_layer);
            sum_fi += fi;
            mxy_I += fi * c_ix[i] * c_iy[i];
        } });

    mxy_I /= sum_fi;

    switch (static_cast<Boundary>(node_in))
    {
    case Boundary::East:
    {
        Constants<Boundary::East> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        sim.mom[momIdx<MomentId::ux>(index)] = 0.f;
        sim.mom[momIdx<MomentId::uy>(index)] = 0.f;
        sim.mom[momIdx<MomentId::mxx>(index)] = 0.f;
        sim.mom[momIdx<MomentId::myy>(index)] = 0.f;
        break;
    }
    case Boundary::North:
    {
        Constants<Boundary::North> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        sim.mom[momIdx<MomentId::ux>(index)] = u_max;
        sim.mom[momIdx<MomentId::uy>(index)] = 0.f;
        sim.mom[momIdx<MomentId::mxx>(index)] = 0.f;
        sim.mom[momIdx<MomentId::myy>(index)] = 0.f;
        break;
    }
    case Boundary::West:
    {
        Constants<Boundary::West> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        sim.mom[momIdx<MomentId::ux>(index)] = 0.f;
        sim.mom[momIdx<MomentId::uy>(index)] = 0.f;
        sim.mom[momIdx<MomentId::mxx>(index)] = 0.f;
        sim.mom[momIdx<MomentId::myy>(index)] = 0.f;
        break;
    }
    case Boundary::South:
    {
        Constants<Boundary::South> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        sim.mom[momIdx<MomentId::ux>(index)] = 0.f;
        sim.mom[momIdx<MomentId::uy>(index)] = 0.f;
        sim.mom[momIdx<MomentId::mxx>(index)] = 0.f;
        sim.mom[momIdx<MomentId::myy>(index)] = 0.f;
        break;
    }
    case Boundary::Northeast:
    {
        Constants<Boundary::Northeast> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        sim.mom[momIdx<MomentId::ux>(index)] = u_max;
        sim.mom[momIdx<MomentId::uy>(index)] = 0.f;
        sim.mom[momIdx<MomentId::mxx>(index)] = 0.f;
        sim.mom[momIdx<MomentId::myy>(index)] = 0.f;
        break;
    }
    case Boundary::Northwest:
    {
        Constants<Boundary::Northwest> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        sim.mom[momIdx<MomentId::ux>(index)] = u_max;
        sim.mom[momIdx<MomentId::uy>(index)] = 0.f;
        sim.mom[momIdx<MomentId::mxx>(index)] = 0.f;
        sim.mom[momIdx<MomentId::myy>(index)] = 0.f;
        break;
    }
    case Boundary::Southwest:
    {
        Constants<Boundary::Southwest> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        sim.mom[momIdx<MomentId::ux>(index)] = 0.f;
        sim.mom[momIdx<MomentId::uy>(index)] = 0.f;
        sim.mom[momIdx<MomentId::mxx>(index)] = 0.f;
        sim.mom[momIdx<MomentId::myy>(index)] = 0.f;
        break;
    }
    case Boundary::Southeast:
    {
        Constants<Boundary::Southeast> c(mask_in);
        sim.mom[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
        sim.mom[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * sim.mom[momIdx<MomentId::mxy>(index)]);
        sim.mom[momIdx<MomentId::ux>(index)] = 0.f;
        sim.mom[momIdx<MomentId::uy>(index)] = 0.f;
        sim.mom[momIdx<MomentId::mxx>(index)] = 0.f;
        sim.mom[momIdx<MomentId::myy>(index)] = 0.f;
        break;
    }
    default:
        break;
    }
}