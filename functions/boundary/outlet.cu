#include "outlet.cuh"

__device__ void outlet_north(CInt *I_s, int x, int y,
                             float *mom_in,
                             float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int i = 0; i < 6; i++)
    {
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();
    int neighbourIndex = from_id(x, y, 3);

    float rho = mom_in[momIdx<MomentId::rho>(neighbourIndex)];
    float ux = mom_in[momIdx<MomentId::ux>(neighbourIndex)];
    float uy = mom_in[momIdx<MomentId::uy>(neighbourIndex)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::mxx>(index)] = (6.f * mxxI * rhoI) / (5.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(-6.f * mxyI * rhoI + rho * ux) / (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = -(-rho - 9.f * myyI * rhoI - 3.f * rho * uy) / (6.f * rho);
}

__device__ void outlet_south(CInt *I_s, int x, int y,
                             float *mom_in,
                             float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int i = 0; i < 6; i++)
    {
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();
    int neighbourIndex = from_id(x, y, 3);

    float rho = mom_in[momIdx<MomentId::rho>(neighbourIndex)];
    float ux = mom_in[momIdx<MomentId::ux>(neighbourIndex)];
    float uy = mom_in[momIdx<MomentId::uy>(neighbourIndex)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::mxx>(index)] = (6.f * mxxI * rhoI) / (5.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(-6.f * mxyI * rhoI - rho * ux) / (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = -(-rho - 9.f * myyI * rhoI - 3.f * rho * uy) / (6.f * rho);
}

__device__ void outlet_east(CInt *I_s, int x, int y,
                            float *mom_in,
                            float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < 6; k++)
    {
        int i = I_s[k];

        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    mxxI /= rhoI;
    mxyI /= rhoI;
    myyI /= rhoI;

    int index = grid_id();
    int neighbourIndex = from_id(x, y, 1);

    float rho = mom_in[momIdx<MomentId::rho>(neighbourIndex)];
    float ux = mom_in[momIdx<MomentId::ux>(neighbourIndex)];
    float uy = mom_in[momIdx<MomentId::uy>(neighbourIndex)];

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::ux>(index)] = ux;
    mom_out[momIdx<MomentId::uy>(index)] = uy;
    mom_out[momIdx<MomentId::mxx>(index)] = -(-rho - 9.f * mxxI * rhoI + 3.f * rho * ux) / (6.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(-6.f * mxyI * rhoI + rho * uy) / (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (6.f * myyI * rhoI) / (5.f * rho);
}

__device__ void outlet_west(CInt *I_s, int x, int y,
                            float *mom_in,
                            float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int i = 0; i < 6; i++)
    {
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();
    int neighbourIndex = from_id(x, y, 3);

    float rho = mom_in[momIdx<MomentId::rho>(neighbourIndex)];
    float ux = mom_in[momIdx<MomentId::ux>(neighbourIndex)];
    float uy = mom_in[momIdx<MomentId::uy>(neighbourIndex)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::mxx>(index)] = -(-rho - 9.f * mxxI * rhoI - 3.f * rho * ux) / (6.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(-6.f * mxyI * rhoI - rho * uy) / (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (6.f * myyI * rhoI) / (5.f * rho);
}

__device__ void outlet_northeast(CInt *I_s, int x, int y,
                                 float *mom_in,
                                 float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < 4; k++)
    {
        int i = I_s[k];

        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();
    int neighbourIndex = from_id(x, y, 3);

    float rho = mom_in[momIdx<MomentId::rho>(neighbourIndex)];
    float ux = mom_in[momIdx<MomentId::ux>(neighbourIndex)];
    float uy = mom_in[momIdx<MomentId::uy>(neighbourIndex)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::ux>(index)] = ux;
    mom_out[momIdx<MomentId::uy>(index)] = uy;
    mom_out[momIdx<MomentId::mxx>(index)] = -(2.f * (-rho - 9.f * mxxI * rhoI + 6.f * mxyI * rhoI + 2.f * rho * ux - rho * uy)) / (9.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(7.f * rho + 18.f * mxxI * rhoI - 132.f * mxyI * rhoI + 18.f * myyI * rhoI + 7.f * rho * ux + 7.f * rho * uy) / (27.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = -(2.f * (-rho + 6.f * mxyI * rhoI - 9.f * myyI * rhoI - rho * ux + 2.f * rho * uy)) / (9.f * rho);
}

__device__ void outlet_northwest(CInt *I_s, int x, int y,
                                 float *mom_in,
                                 float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int i = 0; i < 4; i++)
    {
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();
    int neighbourIndex = from_id(x, y, 3);

    float rho = mom_in[momIdx<MomentId::rho>(neighbourIndex)];
    float ux = mom_in[momIdx<MomentId::ux>(neighbourIndex)];
    float uy = mom_in[momIdx<MomentId::uy>(neighbourIndex)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::mxx>(index)] = (2.f * (rho + 9.f * mxxI * rhoI + 6.f * mxyI * rhoI + 2.f * rho * ux + rho * uy)) / (9.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(-7.f * rho - 18.f * mxxI * rhoI - 132.f * mxyI * rhoI - 18.f * myyI * rhoI + 7.f * rho * ux - 7.f * rho * uy) / (27.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (2.f * (rho + 6.f * mxyI * rhoI + 9.f * myyI * rhoI - rho * ux - 2.f * rho * uy)) / (9.f * rho);
}

__device__ void outlet_southeast(CInt *I_s, int x, int y,
                                 float *mom_in,
                                 float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int k = 0; k < 4; k++)
    {
        int i = I_s[k];

        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();
    int neighbourIndex = from_id(x, y, 3);

    float rho = mom_in[momIdx<MomentId::rho>(neighbourIndex)];
    float ux = mom_in[momIdx<MomentId::ux>(neighbourIndex)];
    float uy = mom_in[momIdx<MomentId::uy>(neighbourIndex)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = rho;
    mom_out[momIdx<MomentId::ux>(index)] = ux;
    mom_out[momIdx<MomentId::uy>(index)] = uy;
    mom_out[momIdx<MomentId::mxx>(index)] = -(2.f * (-rho - 9.f * mxxI * rhoI - 6.f * mxyI * rhoI + 2.f * rho * ux + rho * uy)) / (9.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(-7.f * rho - 18.f * mxxI * rhoI - 132.f * mxyI * rhoI - 18.f * myyI * rhoI - 7.f * rho * ux + 7.f * rho * uy) / (27.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (2.f * (rho + 6.f * mxyI * rhoI + 9.f * myyI * rhoI + rho * ux + 2.f * rho * uy)) / (9.f * rho);
}

__device__ void outlet_southwest(CInt *I_s, int x, int y,
                                 float *mom_in,
                                 float *mom_out)
{
    float rhoI = 0.f;
    float mxxI = 0.f;
    float mxyI = 0.f;
    float myyI = 0.f;

    for (int i = 0; i < 4; i++)
    {
        int indexFrom = from_id(x, y, i);

        float fi = f_i(indexFrom, i, mom_in);

        rhoI += fi;
        mxxI += fi * (c_ix[i] * c_ix[i] - inv_as2);
        mxyI += fi * c_ix[i] * c_iy[i];
        myyI += fi * (c_iy[i] * c_iy[i] - inv_as2);
    }
    int index = grid_id();
    int neighbourIndex = from_id(x, y, 3);

    float rho = mom_in[momIdx<MomentId::rho>(neighbourIndex)];
    float ux = mom_in[momIdx<MomentId::ux>(neighbourIndex)];
    float uy = mom_in[momIdx<MomentId::uy>(neighbourIndex)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::mxx>(index)] = (2.f * (rho + 9.f * mxxI * rhoI - 6.f * mxyI * rhoI + 2.f * rho * ux - rho * uy)) / (9.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = -(7.f * rho + 18.f * mxxI * rhoI - 132.f * mxyI * rhoI + 18.f * myyI * rhoI - 7.f * rho * ux - 7.f * rho * uy) / (27.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = -(2.f * (-rho + 6.f * mxyI * rhoI - 9.f * myyI * rhoI + rho * ux - 2.f * rho * uy)) / (9.f * rho);
}