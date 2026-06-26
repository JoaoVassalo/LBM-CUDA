#include "boundary.cuh"

__device__ void center(moments sim, int x, int y)
{
    float rho = 0.f;
    float ux = 0.f;
    float uy = 0.f;
    float mxx = 0.f;
    float mxy = 0.f;
    float myy = 0.f;

    float f = 0.f;

    for (int i = 1; i < D2Q9::Q; i++)
    {
        f += fi(x, (int)i, sim.layer);

        rho += f;
        ux += f * D2Q9::cx(i);
        uy += f * D2Q9::cy(i);
        mxx += f * (D2Q9::cx(i) * D2Q9::cx(i) - D2Q9::inv_as2);
        mxy += f * (D2Q9::cx(i) * D2Q9::cy(i));
        myy += f * (D2Q9::cy(i) * D2Q9::cy(i) - D2Q9::inv_as2);
    }

    ux /= rho;
    uy /= rho;
    mxx /= rho;
    mxy /= rho;
    myy /= rho;

    sim.layer[1][layerIdx<momId::rho>(x)] = rho;
    sim.layer[1][layerIdx<momId::ux>(x)] = ux;
    sim.layer[1][layerIdx<momId::uy>(x)] = uy;
    sim.layer[1][layerIdx<momId::mxx>(x)] = mxx;
    sim.layer[1][layerIdx<momId::mxy>(x)] = mxy;
    sim.layer[1][layerIdx<momId::myy>(x)] = myy;
}

__device__ void boundary(moments sim, Grid2D grid, int x, int y)
{
    int index = gridId(x, y);

    float rho = 0.f;
    float ux = 0.f;
    float uy = 0.f;
    float mxx = 0.f;
    float mxy = 0.f;
    float myy = 0.f;

    if (y == Geometry::NX - 1)
    {
        ux = physics::u_max;
    }

    float sumfI = 0.f;
    float mxy_I = 0.f;
    float C1 = 0.f;
    float C2 = 0.f;
    float C3 = 0.f;
    float C4 = 0.f;

    float f = 0.f;

    f = fi(x, 0, sim.layer);
    sumfI += f;
    mxy_I = f * (D2Q9::cx(0) * D2Q9::cy(0));
    // O_s
    C1 = D2Q9::w(0) *
             (1.f +
              D2Q9::a_s2 * ux * D2Q9::cx(0) +
              D2Q9::a_s2 * uy * D2Q9::cy(0) +
              D2Q9::a_s4 * 0.5f * ux * ux +
              D2Q9::a_s4 * 0.5f * uy * uy) +
         physics::omega *
             (D2Q9::w(0) * D2Q9::a_s4 * ux * uy *
              (D2Q9::cx(0) * D2Q9::cy(0)));
    C2 = (1.f - physics::omega) *
         (D2Q9::w(0) * D2Q9::a_s4 *
          (D2Q9::cx(0) * D2Q9::cy(0)));

    // I_s
    C3 = D2Q9::w(0) *
         (1.f +
          D2Q9::a_s2 * ux * D2Q9::cx(0) +
          D2Q9::a_s2 * uy * D2Q9::cy(0) +
          D2Q9::a_s4 * 0.5f * ux * ux +
          D2Q9::a_s4 * 0.5f * uy * uy) *
         (D2Q9::cx(0) * D2Q9::cy(0));
    C4 = D2Q9::w(0) *
         D2Q9::a_s4 *
         (D2Q9::cx(0) * D2Q9::cy(0)) *
         (D2Q9::cx(0) * D2Q9::cy(0));

    for (int i = 1; i < D2Q9::Q; i++)
        if (grid.mask[index] & (1u << (i - 1)))
        {
            f = fi(x, i, sim.layer);
            sumfI += f;
            mxy_I += f * (D2Q9::cx(i) * D2Q9::cy(i));

            // I_s
            C3 = D2Q9::w(i) *
                 (1.f +
                  D2Q9::a_s2 * ux * D2Q9::cx(i) +
                  D2Q9::a_s2 * uy * D2Q9::cy(i) +
                  D2Q9::a_s4 * 0.5f * ux * ux +
                  D2Q9::a_s4 * 0.5f * uy * uy) *
                 (D2Q9::cx(i) * D2Q9::cy(i));
            C4 = D2Q9::w(i) *
                 D2Q9::a_s4 *
                 (D2Q9::cx(i) * D2Q9::cy(i)) *
                 (D2Q9::cx(i) * D2Q9::cy(i));

            // O_s
            int k = D2Q9::income(i);

            C1 = D2Q9::w(k) *
                     (1.f +
                      D2Q9::a_s2 * ux * D2Q9::cx(k) +
                      D2Q9::a_s2 * uy * D2Q9::cy(k) +
                      D2Q9::a_s4 * 0.5f * ux * ux +
                      D2Q9::a_s4 * 0.5f * uy * uy) +
                 physics::omega *
                     (D2Q9::w(k) * D2Q9::a_s4 * ux * uy *
                      (D2Q9::cx(k) * D2Q9::cy(k)));
            C2 = (1.f - physics::omega) *
                 (D2Q9::w(k) * D2Q9::a_s4 *
                  (D2Q9::cx(k) * D2Q9::cy(k)));
        }

    mxy = (C3 - C1 * mxy_I) / (C2 * mxy_I - C4);

    rho = sumfI / (C1 + C2 * mxy);

    int id = (y != 0) + (y == Geometry::NY - 1);

    sim.layer[id][layerIdx<momId::rho>(x)] = rho;
    sim.layer[id][layerIdx<momId::ux>(x)] = ux;
    sim.layer[id][layerIdx<momId::uy>(x)] = uy;
    sim.layer[id][layerIdx<momId::mxx>(x)] = mxx;
    sim.layer[id][layerIdx<momId::mxy>(x)] = mxy;
    sim.layer[id][layerIdx<momId::myy>(x)] = myy;
}

__device__ void applyBoundary(moments sim, Grid2D grid, int x, int y)
{
    int index = gridId(x, y);

    if (grid.node[index] != to_u8(Boundary::Center))
        boundary(sim, grid, x, y);
    else
        center(sim, x, y);
}