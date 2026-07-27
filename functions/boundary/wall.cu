#include "../../presets/config.h"
#include "../../presets/physics.h"

#include "../grid_id.cuh"
#include "../from_id.cuh"
#include "../equations/f_i.cuh"

__device__ void wall(CInt *I_s, CInt *O_s, int x, int y,
                     float *mom_in,
                     float *mom_out)
{
    float sum_fi = 0.f;
    float mxy_I = 0.f;

    float Is_up = 0.f;
    float Is_down = 0.f;
    float Os_up = 0.f;
    float Os_down = 0.f;

    int index = grid_id();

#pragma unroll
    for (int k = 0; k < 6; k++)
    {
        int i = I_s[k];
        int index_from = from_id(x, y, i);

        float fi = f_i(index_from, i, mom_in);
        sum_fi += fi;

        mxy_I += fi * c_ix[i] * c_iy[i];

        Is_up += w[i] * c_ix[i] * c_iy[i];

        Is_down += w[i] * a_s4 * c_ix[i] * c_ix[i] * c_iy[i] * c_iy[i];

        /*---------------------------------------------------------------------*/
        i = O_s[k];
        index_from = from_id(x, y, i);

        Os_up += w[i];

        Os_down += w[i] * a_s4 * c_ix[i] * c_iy[i];
    }
    mxy_I /= sum_fi;

    mom_out[momIdx<MomentId::mxy>(index)] = (Is_up - mxy_I * Os_up) / (mxy_I * (1.f - omega) * Os_down - Is_down);              // mxy
    mom_out[momIdx<MomentId::rho>(index)] = sum_fi / ((1.f - omega) * Os_down * mom_out[momIdx<MomentId::mxy>(index)] + Os_up); // rho
}

__device__ void wall_north(int size, CInt *I_s, CInt *O_s,
                           int x, int y,
                           float *mom_in,
                           float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < size; k++)
    {
        //  ENTRADA
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);

        //  SAIDA
        i = O_s[k];

        indexFrom = from_id(x, y, i);

        fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
    }
    int index = grid_id();

    float rho = (3.f * (-4.f * rhoI - 3.f * myyI * rhoI + 3.f * myyI * omega * rhoI)) /
                (-9.f - omega + 3.f * uy + 3.f * omega * uy + 6.f * omega * uy * uy);
    float ux = 0.f;
    float uy = 0.f;

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::mxx>(index)] = (6.f * mxxI * rhoI) / (5.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = (-6.f * mxyI + rho * ux) / (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (-rho - 9.f * myyI * rhoI + 3.f * rho * uy) / (6.f * rho);
}

__device__ void wall_south(int size, CInt *I_s, CInt *O_s,
                           int x, int y,
                           float *mom_in,
                           float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < size; k++)
    {
        //  ENTRADA
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);

        //  SAIDA
        i = O_s[k];

        indexFrom = from_id(x, y, i);

        fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
    }
    int index = grid_id();

    float rho = (3.f * (-4.f * rhoI - 3.f * myyI * rhoI + 3.f * myyI * omega * rhoI)) /
                (-9.f - omega - 3.f * uy - 3.f * omega * uy + 6.f * omega * uy * uy);
    float ux = 0.f;
    float uy = 0.f;

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::mxx>(index)] = (6.f * mxxI * rhoI) / (5.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = (-6.f * mxyI * rhoI - rho * ux) / (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (-rho - 9.f * myyI * rhoI - 3.f * rho * uy) / (6.f * rho);
}

__device__ void wall_east(int size, CInt *I_s, CInt *O_s,
                          int x, int y,
                          float *mom_in,
                          float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < size; k++)
    {
        //  ENTRADA
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);

        //  SAIDA
        i = O_s[k];

        indexFrom = from_id(x, y, i);

        fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
    }
    int index = grid_id();

    float rho = (3.f * (-4.f * rhoI - 3.f * mxxI * rhoI + 3.f * mxxI * omega * rhoI)) /
                (-9.f - omega + 3.f * ux + 3.f * omega * ux + 6.f * omega * ux * ux);
    float ux = 0.f;
    float uy = 0.f;

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::mxx>(index)] =
        (-rho - 9.f * mxxI * rhoI + 3.f * rho * ux) /
        (6.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] =
        (-6.f * mxyI * rhoI + rho * uy) /
        (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] =
        (6.f * myyI * rhoI) /
        (5.f * rho);
}

__device__ void wall_west(int size, CInt *I_s, CInt *O_s,
                          int x, int y,
                          float *mom_in,
                          float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < size; k++)
    {
        //  ENTRADA
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);

        //  SAIDA
        i = O_s[k];

        indexFrom = from_id(x, y, i);

        fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
    }
    int index = grid_id();

    float rho = (3.f * (-4.f * rhoI - 3.f * mxxI * rhoI + 3.f * mxxI * omega * rhoI)) /
                (-9.f - omega - 3.f * ux - 3.f * omega * ux + 6.f * omega * ux * ux);
    float ux = 0.f;
    float uy = 0.f;

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::mxx>(index)] = (-rho - 9.f * mxxI * rhoI - 3.f * rho * ux) /
                                            (6.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = (-6.f * mxyI * rhoI - rho * uy) /
                                            (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (6.f * myyI * rhoI) /
                                            (5.f * rho);
}

__device__ void wall_northeast(int size, CInt *I_s, CInt *O_s,
                               int x, int y,
                               float *mom_in,
                               float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < size; k++)
    {
        //  ENTRADA
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);

        //  SAIDA
        i = O_s[k];

        indexFrom = from_id(x, y, i);

        fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
    }
    int index = grid_id();

    float rho = (12.f * (-3.f * rhoI - 3.f * mxxI * rhoI + 7.f * mxyI * rhoI -
                         3.f * myyI * rhoI + 3.f * mxxI * omega * rhoI -
                         7.f * mxyI * omega * rhoI + 3.f * myyI * omega * rhoI)) /
                (-16.f - 9.f * omega + 14.f * ux + omega * ux +
                 15.f * omega * ux * ux + 14.f * uy + omega * uy -
                 9.f * omega * ux * uy + 15.f * omega * uy * uy);
    float ux = 0.f;
    float uy = 0.f;

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::mxx>(index)] =
        (2.f * (-rho - 9.f * mxxI * rhoI + 6.f * mxyI * rhoI +
                2.f * rho * ux - rho * uy)) /
        (9.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] =
        (7.f * rho + 18.f * mxxI * rhoI - 132.f * mxyI * rhoI +
         18.f * myyI * rhoI + 7.f * rho * ux + 7.f * rho * uy) /
        (27.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] =
        (2.f * (-rho + 6.f * mxyI * rhoI - 9.f * myyI * rhoI -
                rho * ux + 2.f * rho * uy)) /
        (9.f * rho);
}

__device__ void wall_northwest(int size, CInt *I_s, CInt *O_s,
                               int x, int y,
                               float *mom_in,
                               float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < size; k++)
    {
        //  ENTRADA
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);

        //  SAIDA
        i = O_s[k];

        indexFrom = from_id(x, y, i);

        fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
    }
    int index = grid_id();

    float rho = (12.f * (-3.f - 3.f * mxxI - 7.f * mxyI - 3.f * myyI + 3.f * mxxI * omega + 7.f * mxyI * omega + 3.f * myyI * omega) * rhoI) /
                (-16.f - 9.f * omega - 14.f * ux - omega * ux +
                 15.f * omega * ux * ux + 14.f * uy + omega * uy +
                 9.f * omega * ux * uy + 15.f * omega * uy * uy);
    float ux = 0.f;
    float uy = 0.f;

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::mxx>(index)] =
        (2.f * (rho + 9.f * mxxI * rhoI + 6.f * mxyI * rhoI +
                2.f * rho * ux + rho * uy)) /
        (9.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] =
        (-7.f * rho - 18.f * mxxI * rhoI - 132.f * mxyI * rhoI -
         18.f * myyI * rhoI + 7.f * rho * ux - 7.f * rho * uy) /
        (27.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] =
        (2.f * (rho + 6.f * mxyI * rhoI + 9.f * myyI * rhoI -
                rho * ux - 2.f * rho * uy)) /
        (9.f * rho);
}

__device__ void wall_southeast(int size, CInt *I_s, CInt *O_s,
                               int x, int y,
                               float *mom_in,
                               float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < size; k++)
    {
        //  ENTRADA
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);

        //  SAIDA
        i = O_s[k];

        indexFrom = from_id(x, y, i);

        fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
    }
    int index = grid_id();

    float rho = (12.f * (-3.f - 3.f * mxxI - 7.f * mxyI - 3.f * myyI + 3.f * mxxI * omega + 7.f * mxyI * omega + 3.f * myyI * omega) * rhoI) /
                (-16.f - 9.f * omega + 14.f * ux + omega * ux +
                 15.f * omega * ux * ux - 14.f * uy - omega * uy +
                 9.f * omega * ux * uy + 15.f * omega * uy * uy);
    float ux = 0.f;
    float uy = 0.f;

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::mxx>(index)] =
        (2.f * (-rho - 9.f * mxxI * rhoI - 6.f * mxyI * rhoI +
                2.f * rho * ux + rho * uy)) /
        (9.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] =
        (-7.f * rho - 18.f * mxxI * rhoI - 132.f * mxyI * rhoI -
         18.f * myyI * rhoI - 7.f * rho * ux + 7.f * rho * uy) /
        (27.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] =
        (2.f * (rho + 6.f * mxyI * rhoI + 9.f * myyI * rhoI +
                rho * ux + 2.f * rho * uy)) /
        (9.f * rho);
}

__device__ void wall_southwest(int size, CInt *I_s, CInt *O_s,
                               int x, int y,
                               float *mom_in,
                               float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < size; k++)
    {
        //  ENTRADA
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);

        //  SAIDA
        i = O_s[k];

        indexFrom = from_id(x, y, i);

        fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
    }
    int index = grid_id();

    float rho = (12.f * (-3.f - 3.f * mxxI + 7.f * mxyI - 3.f * myyI + 3.f * mxxI * omega - 7.f * mxyI * omega + 3.f * myyI * omega) * rhoI) /
                (-16.f - 9.f * omega - 14.f * ux - omega * ux +
                 15.f * omega * ux * ux - 14.f * uy - omega * uy -
                 9.f * omega * ux * uy + 15.f * omega * uy * uy);
    float ux = 0.f;
    float uy = 0.f;

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::mxx>(index)] =
        (2.f * (rho + 9.f * mxxI * rhoI - 6.f * mxyI * rhoI +
                2.f * rho * ux - rho * uy)) /
        (9.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] =
        (7.f * rho + 18.f * mxxI * rhoI - 132.f * mxyI * rhoI +
         18.f * myyI * rhoI - 7.f * rho * ux - 7.f * rho * uy) /
        (27.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] =
        (2.f * (-rho + 6.f * mxyI * rhoI - 9.f * myyI * rhoI +
                rho * ux - 2.f * rho * uy)) /
        (9.f * rho);
}