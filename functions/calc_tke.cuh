#include "../presets/physics.h"
#include "../presets/geometry.h"
#include <fstream>

inline void calc_tke(std::ofstream &file, float *ux, float *uy, int t)
{
    float tke = 0;
    for (int i = 0; i < (Nx - 1) * (Ny - 1); i++)
    {
        tke += (ux[i] * ux[i] + uy[i] * uy[i]) / (2.f * (float)tf * u_max * u_max);
    }

    file << t << "," << tke << "\n";
}