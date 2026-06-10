#include "../config/physics.h"
#include "../config/geometry.h"
#include "../lbm/build/build_mom.cuh"
#include "indexing.cuh"
#include <fstream>

inline void calc_tke(std::ofstream &file, float *mom, int t)
{
    float tke = 0;
    for (int y = 0; y < Geometry::Ny - 1; y++)
    {
        for (int x = 0; x < Geometry::Nx - 1; x++)
        {
            const float ux = mom[momIdx<MomentId::ux>(x, y)];
            const float uy = mom[momIdx<MomentId::uy>(x, y)];
            tke += (ux * ux + uy * uy) / (2.f * (float)tf * u_max * u_max);
        }
    }

    file << t << "," << tke << "\n";
}
