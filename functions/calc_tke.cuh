#include "../presets/physics.h"
#include "../presets/geometry.h"
#include <fstream>

inline void calc_tke(std::ofstream &file, varUnit *mom, int t)
{
    varUnit tke = 0;
    for (int i = 0; i < (Nx - 1) * (Ny - 1); i++)
    {
        tke += (mom[i + 1] * mom[i + 1]      // ux
                + mom[i + 2] * mom[i + 2]) / // uy
               (2.f * (varUnit)tf * u_max * u_max);
    }

    file << t << "," << tke << "\n";
}