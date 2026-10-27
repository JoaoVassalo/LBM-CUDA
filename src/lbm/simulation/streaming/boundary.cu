#include "boundary.cuh"


__device__ void wallSouthwest(moments sim, Grid2D grid, int x, int y)
{
    int index = gridId(x, y);
    uint8_t mask = grid.mask[index];

    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    varType f = fi(x, 0, sim.layer);

    rhoI += f;

    mxxI += f * (D2Q9::cx(0) * D2Q9::cx(0) - D2Q9::inv_as2);
    mxyI += f * D2Q9::cx(0) * D2Q9::cy(0);
    myyI += f * (D2Q9::cy(0) * D2Q9::cy(0) - D2Q9::inv_as2);


    while (mask)
    {
        int i = __ffs(mask);

        varType f = fi(x, i, sim.layer);

        rhoI += f;

        mxxI += f * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxyI += f * D2Q9::cx(i) * D2Q9::cy(i);
        myyI += f * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);

        mask &= mask - 1;
    }

    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(12) * (-static_cast<varType>(3) - static_cast<varType>(3) * mxxI + static_cast<varType>(7) * mxyI - static_cast<varType>(3) * myyI + static_cast<varType>(3) * mxxI * physics::omega - static_cast<varType>(7) * mxyI * physics::omega + static_cast<varType>(3) * myyI * physics::omega) * rhoI) /
                  (-static_cast<varType>(16) - static_cast<varType>(9) * physics::omega - static_cast<varType>(14) * ux - physics::omega * ux + static_cast<varType>(15) * physics::omega * ux * ux - static_cast<varType>(14) * uy - physics::omega * uy - static_cast<varType>(9) * physics::omega * ux * uy + static_cast<varType>(15) * physics::omega * uy * uy);

    sim.mom[momIdx<momId::rho>(index)] = rho;
    sim.mom[momIdx<momId::ux>(index)] = ux;
    sim.mom[momIdx<momId::uy>(index)] = uy;
    sim.mom[momIdx<momId::mxx>(index)] = (static_cast<varType>(2) * (rho + static_cast<varType>(9) * mxxI * rhoI - static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(2) * rho * ux - rho * uy)) /
                                         (static_cast<varType>(9) * rho);
    sim.mom[momIdx<momId::mxy>(index)] = -(static_cast<varType>(7) * rho + static_cast<varType>(18) * mxxI * rhoI - static_cast<varType>(132) * mxyI * rhoI + static_cast<varType>(18) * myyI * rhoI - static_cast<varType>(7) * rho * ux - static_cast<varType>(7) * rho * uy) /
                                         (static_cast<varType>(27) * rho);
    sim.mom[momIdx<momId::myy>(index)] = -(static_cast<varType>(2) * (-rho + static_cast<varType>(6) * mxyI * rhoI - static_cast<varType>(9) * myyI * rhoI + rho * ux - static_cast<varType>(2) * rho * uy)) /
                                         (static_cast<varType>(9) * rho);
}

__device__ void wallSoutheast(moments sim, Grid2D grid, int x, int y)
{
    int index = gridId(x, y);
    uint8_t mask = grid.mask[index];

    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    varType f = fi(x, 0, sim.layer);

    rhoI += f;

    mxxI += f * (D2Q9::cx(0) * D2Q9::cx(0) - D2Q9::inv_as2);
    mxyI += f * D2Q9::cx(0) * D2Q9::cy(0);
    myyI += f * (D2Q9::cy(0) * D2Q9::cy(0) - D2Q9::inv_as2);


    while (mask)
    {
        int i = __ffs(mask);

        varType f = fi(x, i, sim.layer);

        rhoI += f;

        mxxI += f * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxyI += f * D2Q9::cx(i) * D2Q9::cy(i);
        myyI += f * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);

        mask &= mask - 1;
    }

    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(12) * (-static_cast<varType>(3) - static_cast<varType>(3) * mxxI - static_cast<varType>(7) * mxyI - static_cast<varType>(3) * myyI + static_cast<varType>(3) * mxxI * physics::omega + static_cast<varType>(7) * mxyI * physics::omega + static_cast<varType>(3) * myyI * physics::omega) * rhoI) /
                  (-static_cast<varType>(16) - static_cast<varType>(9) * physics::omega + static_cast<varType>(14) * ux + physics::omega * ux + static_cast<varType>(15) * physics::omega * ux * ux - static_cast<varType>(14) * uy - physics::omega * uy + static_cast<varType>(9) * physics::omega * ux * uy + static_cast<varType>(15) * physics::omega * uy * uy);

    sim.mom[momIdx<momId::rho>(index)] = rho;
    sim.mom[momIdx<momId::ux>(index)] = ux;
    sim.mom[momIdx<momId::uy>(index)] = uy;
    sim.mom[momIdx<momId::mxx>(index)] = -(static_cast<varType>(2) * (-rho - static_cast<varType>(9) * mxxI * rhoI - static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(2) * rho * ux + rho * uy)) /
                                         (static_cast<varType>(9) * rho);
    sim.mom[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(7) * rho - static_cast<varType>(18) * mxxI * rhoI - static_cast<varType>(132) * mxyI * rhoI - static_cast<varType>(18) * myyI * rhoI - static_cast<varType>(7) * rho * ux + static_cast<varType>(7) * rho * uy) /
                                         (static_cast<varType>(27) * rho);
    sim.mom[momIdx<momId::myy>(index)] = (static_cast<varType>(2) * (rho + static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(9) * myyI * rhoI + rho * ux + static_cast<varType>(2) * rho * uy)) /
                                         (static_cast<varType>(9) * rho);
}

__device__ void wallNorthwest(moments sim, Grid2D grid, int x, int y)
{
    int index = gridId(x, y);
    uint8_t mask = grid.mask[index];

    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    varType f = fi(x, 0, sim.layer);

    rhoI += f;

    mxxI += f * (D2Q9::cx(0) * D2Q9::cx(0) - D2Q9::inv_as2);
    mxyI += f * D2Q9::cx(0) * D2Q9::cy(0);
    myyI += f * (D2Q9::cy(0) * D2Q9::cy(0) - D2Q9::inv_as2);


    while (mask)
    {
        int i = __ffs(mask);

        varType f = fi(x, i, sim.layer);

        rhoI += f;

        mxxI += f * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxyI += f * D2Q9::cx(i) * D2Q9::cy(i);
        myyI += f * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);

        mask &= mask - 1;
    }

    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(12) * (-static_cast<varType>(3) - static_cast<varType>(3) * mxxI - static_cast<varType>(7) * mxyI - static_cast<varType>(3) * myyI + static_cast<varType>(3) * mxxI * physics::omega + static_cast<varType>(7) * mxyI * physics::omega + static_cast<varType>(3) * myyI * physics::omega) * rhoI) /
                  (-static_cast<varType>(16) - static_cast<varType>(9) * physics::omega - static_cast<varType>(14) * ux - physics::omega * ux + static_cast<varType>(15) * physics::omega * ux * ux + static_cast<varType>(14) * uy + physics::omega * uy + static_cast<varType>(9) * physics::omega * ux * uy + static_cast<varType>(15) * physics::omega * uy * uy);

    sim.mom[momIdx<momId::rho>(index)] = rho;
    sim.mom[momIdx<momId::ux>(index)] = ux;
    sim.mom[momIdx<momId::uy>(index)] = uy;
    sim.mom[momIdx<momId::mxx>(index)] = (static_cast<varType>(2) * (rho + static_cast<varType>(9) * mxxI * rhoI + static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(2) * rho * ux + rho * uy)) /
                                         (static_cast<varType>(9) * rho);
    sim.mom[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(7) * rho - static_cast<varType>(18) * mxxI * rhoI - static_cast<varType>(132) * mxyI * rhoI - static_cast<varType>(18) * myyI * rhoI + static_cast<varType>(7) * rho * ux - static_cast<varType>(7) * rho * uy) /
                                         (static_cast<varType>(27) * rho);
    sim.mom[momIdx<momId::myy>(index)] = (static_cast<varType>(2) * (rho + static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(9) * myyI * rhoI - rho * ux - static_cast<varType>(2) * rho * uy)) /
                                         (static_cast<varType>(9) * rho);
}

__device__ void wallNortheast(moments sim, Grid2D grid, int x, int y)
{
    int index = gridId(x, y);
    uint8_t mask = grid.mask[index];

    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    varType f = fi(x, 0, sim.layer);

    rhoI += f;

    mxxI += f * (D2Q9::cx(0) * D2Q9::cx(0) - D2Q9::inv_as2);
    mxyI += f * D2Q9::cx(0) * D2Q9::cy(0);
    myyI += f * (D2Q9::cy(0) * D2Q9::cy(0) - D2Q9::inv_as2);


    while (mask)
    {
        int i = __ffs(mask);

        varType f = fi(x, i, sim.layer);

        rhoI += f;

        mxxI += f * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxyI += f * D2Q9::cx(i) * D2Q9::cy(i);
        myyI += f * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);

        mask &= mask - 1;
    }

    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(12) * (-static_cast<varType>(3) * rhoI - static_cast<varType>(3) * mxxI * rhoI + static_cast<varType>(7) * mxyI * rhoI - static_cast<varType>(3) * myyI * rhoI + static_cast<varType>(3) * mxxI * physics::omega * rhoI - static_cast<varType>(7) * mxyI * physics::omega * rhoI + static_cast<varType>(3) * myyI * physics::omega * rhoI)) /
                  (-static_cast<varType>(16) - static_cast<varType>(9) * physics::omega + static_cast<varType>(14) * ux + physics::omega * ux + static_cast<varType>(15) * physics::omega * ux * ux + static_cast<varType>(14) * uy + physics::omega * uy - static_cast<varType>(9) * physics::omega * ux * uy + static_cast<varType>(15) * physics::omega * uy * uy);

    sim.mom[momIdx<momId::rho>(index)] = rho;
    sim.mom[momIdx<momId::ux>(index)] = ux;
    sim.mom[momIdx<momId::uy>(index)] = uy;
    sim.mom[momIdx<momId::mxx>(index)] = -(static_cast<varType>(2) * (-rho - static_cast<varType>(9) * mxxI * rhoI + static_cast<varType>(6) * mxyI * rhoI + static_cast<varType>(2) * rho * ux - rho * uy)) /
                                         (static_cast<varType>(9) * rho);
    sim.mom[momIdx<momId::mxy>(index)] = -(static_cast<varType>(7) * rho + static_cast<varType>(18) * mxxI * rhoI - static_cast<varType>(132) * mxyI * rhoI + static_cast<varType>(18) * myyI * rhoI + static_cast<varType>(7) * rho * ux + static_cast<varType>(7) * rho * uy) /
                                         (static_cast<varType>(27) * rho);
    sim.mom[momIdx<momId::myy>(index)] = -(static_cast<varType>(2) * (-rho + static_cast<varType>(6) * mxyI * rhoI - static_cast<varType>(9) * myyI * rhoI - rho * ux + static_cast<varType>(2) * rho * uy)) /
                                         (static_cast<varType>(9) * rho);
}

__device__ void wallSouth(moments sim, Grid2D grid, int x, int y)
{
    int index = gridId(x, y);
    uint8_t mask = grid.mask[index];

    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    varType f = fi(x, 0, sim.layer);

    rhoI += f;

    mxxI += f * (D2Q9::cx(0) * D2Q9::cx(0) - D2Q9::inv_as2);
    mxyI += f * D2Q9::cx(0) * D2Q9::cy(0);
    myyI += f * (D2Q9::cy(0) * D2Q9::cy(0) - D2Q9::inv_as2);


    while (mask)
    {
        int i = __ffs(mask);

        varType f = fi(x, i, sim.layer);

        rhoI += f;

        mxxI += f * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxyI += f * D2Q9::cx(i) * D2Q9::cy(i);
        myyI += f * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);

        mask &= mask - 1;
    }

    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(3) * (-static_cast<varType>(4) * rhoI - static_cast<varType>(3) * myyI * rhoI + static_cast<varType>(3) * myyI * physics::omega * rhoI)) /
                  (-static_cast<varType>(9) - physics::omega - static_cast<varType>(3) * uy - static_cast<varType>(3) * physics::omega * uy + static_cast<varType>(6) * physics::omega * uy * uy);

    sim.mom[momIdx<momId::rho>(index)] = rho;
    sim.mom[momIdx<momId::ux>(index)] = ux;
    sim.mom[momIdx<momId::uy>(index)] = uy;
    sim.mom[momIdx<momId::mxx>(index)] = (static_cast<varType>(6) * mxxI * rhoI) / (static_cast<varType>(5) * rho);
    sim.mom[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(6) * mxyI * rhoI - rho * ux) / (static_cast<varType>(3) * rho);
    sim.mom[momIdx<momId::myy>(index)] = -(-rho - static_cast<varType>(9) * myyI * rhoI - static_cast<varType>(3) * rho * uy) / (static_cast<varType>(6) * rho);
}

__device__ void wallNorth(moments sim, Grid2D grid, int x, int y)
{
    int index = gridId(x, y);
    uint8_t mask = grid.mask[index];

    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    varType f = fi(x, 0, sim.layer);

    rhoI += f;

    mxxI += f * (D2Q9::cx(0) * D2Q9::cx(0) - D2Q9::inv_as2);
    mxyI += f * D2Q9::cx(0) * D2Q9::cy(0);
    myyI += f * (D2Q9::cy(0) * D2Q9::cy(0) - D2Q9::inv_as2);


    while (mask)
    {
        int i = __ffs(mask);

        varType f = fi(x, i, sim.layer);

        rhoI += f;

        mxxI += f * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxyI += f * D2Q9::cx(i) * D2Q9::cy(i);
        myyI += f * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);

        mask &= mask - 1;
    }

    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    varType ux = static_cast<varType>(physics::u_max);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(3) * (-static_cast<varType>(4) * rhoI - static_cast<varType>(3) * myyI * rhoI + static_cast<varType>(3) * myyI * physics::omega * rhoI)) /
                  (-static_cast<varType>(9) - physics::omega + static_cast<varType>(3) * uy + static_cast<varType>(3) * physics::omega * uy + static_cast<varType>(6) * physics::omega * uy * uy);

    sim.mom[momIdx<momId::rho>(index)] = rho;
    sim.mom[momIdx<momId::ux>(index)] = ux;
    sim.mom[momIdx<momId::uy>(index)] = uy;
    sim.mom[momIdx<momId::mxx>(index)] = (static_cast<varType>(6) * mxxI * rhoI) / (static_cast<varType>(5) * rho);
    sim.mom[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(6) * mxyI * rhoI + rho * ux) / (static_cast<varType>(3) * rho);
    sim.mom[momIdx<momId::myy>(index)] = -(-rho - static_cast<varType>(9) * myyI * rhoI + static_cast<varType>(3) * rho * uy) / (static_cast<varType>(6) * rho);
}

__device__ void wallWest(moments sim, Grid2D grid, int x, int y)
{
    int index = gridId(x, y);
    uint8_t mask = grid.mask[index];

    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    varType f = fi(x, 0, sim.layer);

    rhoI += f;

    mxxI += f * (D2Q9::cx(0) * D2Q9::cx(0) - D2Q9::inv_as2);
    mxyI += f * D2Q9::cx(0) * D2Q9::cy(0);
    myyI += f * (D2Q9::cy(0) * D2Q9::cy(0) - D2Q9::inv_as2);


    while (mask)
    {
        int i = __ffs(mask);

        varType f = fi(x, i, sim.layer);

        rhoI += f;

        mxxI += f * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxyI += f * D2Q9::cx(i) * D2Q9::cy(i);
        myyI += f * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);

        mask &= mask - 1;
    }

    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(3) * (-static_cast<varType>(4) * rhoI - static_cast<varType>(3) * mxxI * rhoI + static_cast<varType>(3) * mxxI * physics::omega * rhoI)) /
                  (-static_cast<varType>(9) - physics::omega - static_cast<varType>(3) * ux - static_cast<varType>(3) * physics::omega * ux + static_cast<varType>(6) * physics::omega * ux * ux);

    sim.mom[momIdx<momId::rho>(index)] = rho;
    sim.mom[momIdx<momId::ux>(index)] = ux;
    sim.mom[momIdx<momId::uy>(index)] = uy;
    sim.mom[momIdx<momId::mxx>(index)] = -(-rho - static_cast<varType>(9) * mxxI * rhoI - static_cast<varType>(3) * rho * ux) /
                                         (static_cast<varType>(6) * rho);
    sim.mom[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(6) * mxyI * rhoI - rho * uy) /
                                         (static_cast<varType>(3) * rho);
    sim.mom[momIdx<momId::myy>(index)] = (static_cast<varType>(6) * myyI * rhoI) /
                                         (static_cast<varType>(5) * rho);
}

__device__ void wallEast(moments sim, Grid2D grid, int x, int y)
{
    int index = gridId(x, y);
    uint8_t mask = grid.mask[index];

    varType rhoI = static_cast<varType>(0);
    varType mxxI = static_cast<varType>(0);
    varType mxyI = static_cast<varType>(0);
    varType myyI = static_cast<varType>(0);

    varType f = fi(x, 0, sim.layer);

    rhoI += f;

    mxxI += f * (D2Q9::cx(0) * D2Q9::cx(0) - D2Q9::inv_as2);
    mxyI += f * D2Q9::cx(0) * D2Q9::cy(0);
    myyI += f * (D2Q9::cy(0) * D2Q9::cy(0) - D2Q9::inv_as2);


    while (mask)
    {
        int i = __ffs(mask);

        varType f = fi(x, i, sim.layer);

        rhoI += f;

        mxxI += f * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxyI += f * D2Q9::cx(i) * D2Q9::cy(i);
        myyI += f * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);

        mask &= mask - 1;
    }

    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    varType ux = static_cast<varType>(0);
    varType uy = static_cast<varType>(0);

    varType rho = (static_cast<varType>(3) * (-static_cast<varType>(4) * rhoI - static_cast<varType>(3) * mxxI * rhoI + static_cast<varType>(3) * mxxI * physics::omega * rhoI)) /
                  (-static_cast<varType>(9) - physics::omega + static_cast<varType>(3) * ux + static_cast<varType>(3) * physics::omega * ux + static_cast<varType>(6) * physics::omega * ux * ux);

    sim.mom[momIdx<momId::rho>(index)] = rho;
    sim.mom[momIdx<momId::ux>(index)] = ux;
    sim.mom[momIdx<momId::uy>(index)] = uy;
    sim.mom[momIdx<momId::mxx>(index)] = -(-rho - static_cast<varType>(9) * mxxI * rhoI + static_cast<varType>(3) * rho * ux) / (static_cast<varType>(6) * rho);
    sim.mom[momIdx<momId::mxy>(index)] = -(-static_cast<varType>(6) * mxyI * rhoI + rho * uy) / (static_cast<varType>(3) * rho);
    sim.mom[momIdx<momId::myy>(index)] = (static_cast<varType>(6) * myyI * rhoI) / (static_cast<varType>(5) * rho);
}
__device__ void center(moments sim, int x, int y)
{
    int index = gridId(x, y);
    varType rho = 0.;
    varType ux = 0.;
    varType uy = 0.;
    varType mxx = 0.;
    varType mxy = 0.;
    varType myy = 0.;

    constexprFor<0, D2Q9::Q>([&](auto i)
    {
        constexpr int direction = decltype(i)::value;
        const varType f = fi(x, direction, sim.layer);
        rho += f;
        ux += f * D2Q9::cx(direction);
        uy += f * D2Q9::cy(direction);
        mxx += f * (D2Q9::cx(direction) * D2Q9::cx(direction) - D2Q9::inv_as2);
        mxy += f * (D2Q9::cx(direction) * D2Q9::cy(direction));
        myy += f * (D2Q9::cy(direction) * D2Q9::cy(direction) - D2Q9::inv_as2);
    });

    ux /= rho;
    uy /= rho;
    mxx /= rho;
    mxy /= rho;
    myy /= rho;

    sim.mom[momIdx<momId::rho>(index)] = rho;
    sim.mom[momIdx<momId::ux>(index)] = ux;
    sim.mom[momIdx<momId::uy>(index)] = uy;
    sim.mom[momIdx<momId::mxx>(index)] = mxx;
    sim.mom[momIdx<momId::mxy>(index)] = mxy;
    sim.mom[momIdx<momId::myy>(index)] = myy;
}

__device__ void applyBoundary(moments sim, Grid2D grid, int x, int y)
{
    if (x == 0 && y == 0)
        wallSouthwest(sim, grid, x, y);
    else if (x == Geometry::NX - 1 && y == 0)
        wallSoutheast(sim, grid, x, y);
    else if (x == 0 && y == Geometry::NY - 1)
        wallNorthwest(sim, grid, x, y);
    else if (x == Geometry::NX - 1 && y == Geometry::NY - 1)
        wallNortheast(sim, grid, x, y);
    else if (y == 0)
        wallSouth(sim, grid, x, y);
    else if (y == Geometry::NY - 1)
        wallNorth(sim, grid, x, y);
    else if (x == 0)
        wallWest(sim, grid, x, y);
    else if (x == Geometry::NX - 1)
        wallEast(sim, grid, x, y);
    else
        center(sim, x, y);
}
