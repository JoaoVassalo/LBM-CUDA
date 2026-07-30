#include "boundary.cuh"

__device__ void boundary(int size, CInt *Is, CInt *Os, int x, int y, varUnit *mom_in,
                         varUnit *mom_out)
{
    int index = grid_id();

    varUnit rho = 0.f;
    varUnit ux = 0.f;
    varUnit uy = 0.f;
    varUnit mxx = 0.f;
    varUnit mxy = 0.f;
    varUnit myy = 0.f;

    varUnit sumfI = 0.f;
    varUnit mxy_I = 0.f;
    varUnit C1 = 0.f;
    varUnit C2 = 0.f;
    varUnit C3 = 0.f;
    varUnit C4 = 0.f;

    varUnit f = 0.f;

    for (int k = 0; k < size; k++)
    {
        int i = Is[k];
        f = f_i(index, i, mom_in);
        sumfI += f;
        mxy_I += f * (c_ix[i] * c_iy[i]);

        C3 = w[i] *
             (1.f +
              a_s2 * ux * c_ix[i] +
              a_s2 * uy * c_iy[i] +
              a_s4 * 0.5f * ux * ux * (c_ix[i] * c_ix[i] - inv_as2) +
              a_s4 * 0.5f * uy * uy * (c_iy[i] * c_iy[i] - inv_as2)) *
             (c_ix[i] * c_iy[i]);
        C4 = w[i] *
             a_s4 *
             (c_ix[i] * c_iy[i]) *
             (c_ix[i] * c_iy[i]);

        // O_s
        i = Os[k];

        C1 = w[i] *
                 (1.f +
                  a_s2 * ux * c_ix[i] +
                  a_s2 * uy * c_iy[i] +
                  a_s4 * 0.5f * ux * ux * (c_ix[i] * c_ix[i] - inv_as2) +
                  a_s4 * 0.5f * uy * uy * (c_iy[i] * c_iy[i] - inv_as2)) +
             omega *
                 (w[i] * a_s4 * ux * uy *
                  (c_ix[i] * c_iy[i]));
        C2 = (1.f - omega) *
             (w[i] * a_s4 *
              (c_ix[i] * c_iy[i]));
    }

    mxy = (C3 - C1 * mxy_I) / (C2 * mxy_I - C4);

    rho = sumfI / (C1 + C2 * mxy);

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::ux>(index)] = ux;
    mom_out[momIdx<MomentId::uy>(index)] = uy;
    mom_out[momIdx<MomentId::mxx>(index)] = mxx;
    mom_out[momIdx<MomentId::mxy>(index)] = mxy;
    mom_out[momIdx<MomentId::myy>(index)] = myy;
}