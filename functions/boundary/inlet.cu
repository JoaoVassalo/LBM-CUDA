#include "inlet.cuh"

#include <cstdio>

__device__ void inlet_north(CInt *I_s, int x, int y,
                            varUnit *mom_in,
                            varUnit *mom_out)
{
    varUnit rhoI = 0.f;
    varUnit mxxI = 0.f;
    varUnit mxyI = 0.f;
    varUnit myyI = 0.f;

    for (int i = 0; i < 6; i++)
    {
        int indexFrom = from_id(x, y, i);

        varUnit fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();

    varUnit rho = mom_in[momIdx<MomentId::rho>(index)];
    varUnit ux = u_max;
    varUnit uy = 0.f;
    varUnit mxx = mom_in[momIdx<MomentId::mxx>(index)];
    varUnit mxy = mom_in[momIdx<MomentId::mxy>(index)];
    varUnit myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = ((4.f + 3.f * myyI) * rhoI) / 3.f * (1.f + uy);
    mom_out[momIdx<MomentId::mxx>(index)] = (18.f * mxxI * (1.f + uy)) / (20.f + 15.f * myyI);
    mom_out[momIdx<MomentId::mxy>(index)] = (-18.f * mxyI * (1.f + uy) + 4.f * ux + 3.f * myyI * ux) / (-12.f - 9.f * myyI);
    mom_out[momIdx<MomentId::myy>(index)] = (2.f + 15.f * myyI - 6.f * uy + 9.f * myyI * uy) / (12.f + 9.f * myyI);
}

__device__ void inlet_south(CInt *I_s, int x, int y,
                            varUnit *mom_in,
                            varUnit *mom_out)
{
    varUnit rhoI = 0.f;
    varUnit mxxI = 0.f;
    varUnit mxyI = 0.f;
    varUnit myyI = 0.f;

    for (int i = 0; i < 6; i++)
    {
        int indexFrom = from_id(x, y, i);

        varUnit fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();

    varUnit rho = mom_in[momIdx<MomentId::rho>(index)];
    varUnit ux = u_max;
    varUnit uy = 0.f;
    varUnit mxx = mom_in[momIdx<MomentId::mxx>(index)];
    varUnit mxy = mom_in[momIdx<MomentId::mxy>(index)];
    varUnit myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = ((-4.f - 3.f * myyI) * rhoI) / 3.f * (-1.f + uy);
    mom_out[momIdx<MomentId::mxx>(index)] = (18.f * mxxI * (-1.f + uy)) / (20.f + 15.f * myyI);
    mom_out[momIdx<MomentId::mxy>(index)] = (18.f * mxyI * (1.f - uy) + 4.f * ux + 3.f * myyI * ux) / (12.f + 9.f * myyI);
    mom_out[momIdx<MomentId::myy>(index)] = (2.f + 15.f * myyI + 6.f * uy - 9.f * myyI * uy) / (12.f + 9.f * myyI);
}

__device__ void inlet_east(CInt *I_s, int x, int y,
                           varUnit *mom_in,
                           varUnit *mom_out)
{
    varUnit rhoI = 0.f;
    varUnit mxxI = 0.f;
    varUnit mxyI = 0.f;
    varUnit myyI = 0.f;

    for (int i = 0; i < 6; i++)
    {
        int indexFrom = from_id(x, y, i);

        varUnit fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();

    varUnit rho = mom_in[momIdx<MomentId::rho>(index)];
    varUnit ux = u_max;
    varUnit uy = 0.f;
    varUnit mxx = mom_in[momIdx<MomentId::mxx>(index)];
    varUnit mxy = mom_in[momIdx<MomentId::mxy>(index)];
    varUnit myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = ((4.f + 3.f * mxxI) * rhoI) / 3.f * (1.f + ux);
    mom_out[momIdx<MomentId::mxx>(index)] = (2.f + 15.f * mxxI - 6.f * ux + 9.f * mxxI * ux) / (12.f + 9.f * mxxI);
    mom_out[momIdx<MomentId::mxy>(index)] = (18.f * mxyI * (1.f + ux) - 4.f * uy - 3.f * mxxI * uy) / (12.f + 9.f * mxxI);
    mom_out[momIdx<MomentId::myy>(index)] = (18.f * myyI * (1.f + ux)) / (20.f + 15.f * myyI);
}

__device__ void inlet_west(CInt *I_s, int x, int y,
                           varUnit *mom_in,
                           varUnit *mom_out)
{
    varUnit rhoI = 0.f;
    varUnit mxxI = 0.f;
    varUnit mxyI = 0.f;
    varUnit myyI = 0.f;

    for (int k = 0; k < 6; k++)
    {
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        varUnit fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    int index = grid_id();

    varUnit ux = u_max;
    varUnit uy = 0.f;

    varUnit rho = (-4.f * rhoI - 3.f * mxxI * rhoI) / (3.f * (-1.f + ux));

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::ux>(index)] = ux;
    mom_out[momIdx<MomentId::uy>(index)] = uy;
    mom_out[momIdx<MomentId::mxx>(index)] = -(-rho - 9.f * mxxI * rhoI - 3.f * rho * ux) / (6.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(-6.f * mxyI * rhoI - rho * uy) / (6.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = -(6 * myyI * rhoI) / (5.f * rho);
}

__device__ void inlet_northeast(CInt *I_s, int x, int y,
                                varUnit *mom_in,
                                varUnit *mom_out)
{
    varUnit rhoI = 0.f;
    varUnit mxxI = 0.f;
    varUnit mxyI = 0.f;
    varUnit myyI = 0.f;

    for (int i = 0; i < 4; i++)
    {
        int indexFrom = from_id(x, y, i);

        varUnit fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();

    varUnit rho = mom_in[momIdx<MomentId::rho>(index)];
    varUnit ux = u_max;
    varUnit uy = 0.f;
    varUnit mxx = mom_in[momIdx<MomentId::mxx>(index)];
    varUnit mxy = mom_in[momIdx<MomentId::mxy>(index)];
    varUnit myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = (3.f * (3.f * rhoI + 3.f * mxxI * rhoI - 7.f * mxyI * rhoI + 3.f * myyI * rhoI)) / 4.f * (1.f + ux + uy);
    mom_out[momIdx<MomentId::mxx>(index)] = (2.f * (1.f + 5.f * mxxI - 5.f * mxyI + myyI - 2.f * ux + 2.f * mxxI * ux + 2.f * mxyI * ux - 2.f * myyI * ux + uy + 5.f * mxxI * uy - 5.f * mxyI * uy + myyI * uy)) / (3 * (3.f + 3.f * mxxI - 7.f * mxyI + 3.f * myyI));
    mom_out[momIdx<MomentId::mxy>(index)] = -(7.f + 15.f * mxxI - 75.f * mxyI + 15.f * myyI) * (1.f + ux + uy) / (9.f * (3.f + 3.f * mxxI - 7.f * mxyI + 3.f * myyI));
    mom_out[momIdx<MomentId::myy>(index)] = (2.f * (1.f + mxxI - 5.f * mxyI + 5.f * myyI + ux + mxxI * ux - 5.f * mxyI * ux + 5.f * myyI * ux - 2.f * uy - 2.f * mxxI * uy + 2.f * mxyI * uy + 2.f * myyI * uy)) / (3.f * (3.f + 3.f * mxxI - 7.f * mxyI + 3.f * myyI));
}

__device__ void inlet_northwest(CInt *I_s, int x, int y,
                                varUnit *mom_in,
                                varUnit *mom_out)
{
    varUnit rhoI = 0.f;
    varUnit mxxI = 0.f;
    varUnit mxyI = 0.f;
    varUnit myyI = 0.f;

    for (int k = 0; k < 4; k++)
    {
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        varUnit fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();

    varUnit ux = u_max;
    varUnit uy = 0.f;

    varUnit rho = -(3.f * (3.f * rhoI + 3.f * mxxI * rhoI + 7.f * mxyI * rhoI + 3.f * myyI * rhoI)) /
                4.f * (-1.f + ux - uy);

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::ux>(index)] = u_max; // ux
    mom_out[momIdx<MomentId::uy>(index)] = 0.f;   // uy
    mom_out[momIdx<MomentId::mxx>(index)] = (2.f * (rho + 9.f * mxxI * rhoI + 6.f * mxyI * rhoI + 2.f * rho * ux + rho * uy)) /
                                            (9.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(-7.f * rho - 18.f * mxxI * rhoI - 132.f * mxyI * rhoI - 18.f * myyI * rhoI + 7.f * rho * ux - 7.f * rho * uy) /
                                            (27.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (2.f * (rho + 6.f * mxyI * rhoI + 9.f * myyI * rhoI - rho * ux - 2.f * rho * uy)) /
                                            (9.f * rho);
}

__device__ void inlet_southeast(CInt *I_s, int x, int y,
                                varUnit *mom_in,
                                varUnit *mom_out)
{
    varUnit rhoI = 0.f;
    varUnit mxxI = 0.f;
    varUnit mxyI = 0.f;
    varUnit myyI = 0.f;

    for (int i = 0; i < 4; i++)
    {
        int indexFrom = from_id(x, y, i);

        varUnit fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();

    varUnit rho = mom_in[momIdx<MomentId::rho>(index)];
    varUnit ux = u_max;
    varUnit uy = 0.f;
    varUnit mxx = mom_in[momIdx<MomentId::mxx>(index)];
    varUnit mxy = mom_in[momIdx<MomentId::mxy>(index)];
    varUnit myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = (3.f * (3.f * rhoI + 3.f * mxxI * rhoI + 7.f * mxyI * rhoI + 3.f * myyI * rhoI)) / (4.f * (1.f + ux - uy));
    mom_out[momIdx<MomentId::mxx>(index)] = (2.f * (1.f + 5.f * mxxI + 5.f * mxyI + myyI - 2.f * ux + 2.f * mxxI * ux - 2.f * mxyI * ux - 2.f * myyI * ux - uy - 5.f * mxxI * uy - 5.f * mxyI * uy - myyI * uy)) / (3.f * (3.f + 3.f * mxxI + 7.f * mxyI + 3.f * myyI));
    mom_out[momIdx<MomentId::mxy>(index)] = ((7.f + 15.f * mxxI + 75.f * mxyI + 15.f * myyI) * (1.f + ux - uy)) / (27.f + 27.f * mxxI + 63.f * mxyI + 27.f * myyI);
    mom_out[momIdx<MomentId::myy>(index)] = (2.f * (1.f + mxxI + 5.f * mxyI + 5.f * myyI + ux + mxxI * ux + 5.f * mxyI * ux + 5.f * myyI * ux + 2.f * uy + 2.f * mxxI * uy + 2.f * mxyI * uy - 2.f * myyI * uy)) / (3.f * (3.f + 3.f * mxxI + 7.f * mxyI + 3.f * myyI));
}

__device__ void inlet_southwest(CInt *I_s, int x, int y,
                                varUnit *mom_in,
                                varUnit *mom_out)
{
    varUnit rhoI = 0.f;
    varUnit mxxI = 0.f;
    varUnit mxyI = 0.f;
    varUnit myyI = 0.f;

    for (int k = 0; k < 4; k++)
    {
        int i = I_s[k];
        int indexFrom = from_id(x, y, i);

        varUnit fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }

    int index = grid_id();

    varUnit ux = u_max;
    varUnit uy = 0.f;

    varUnit rho = -(3.f * (3.f * rhoI + 3.f * mxxI * rhoI - 7.f * mxyI * rhoI + 3.f * myyI * rhoI)) /
                (4.f * (-1.f + ux + uy));
    ;

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::ux>(index)] = ux;
    mom_out[momIdx<MomentId::uy>(index)] = uy;
    mom_out[momIdx<MomentId::mxx>(index)] = (2.f * (rho + 9.f * mxxI * rhoI - 6.f * mxyI * rhoI + 2.f * rho * ux - rho * uy)) /
                                            (9.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(7.f * rho + 18.f * mxxI * rhoI - 132.f * mxyI * rhoI + 18.f * myyI * rhoI - 7.f * rho * ux - 7.f * rho * uy) /
                                            (27.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = -(2.f * (-rho + 6.f * mxyI * rhoI - 9.f * myyI * rhoI + rho * ux - 2.f * rho * uy)) /
                                            (9.f * rho);
}