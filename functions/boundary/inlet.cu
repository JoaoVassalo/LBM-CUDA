#include "inlet.cuh"

__device__ void inlet_north(CInt *I_s, int x, int y,
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

    float rho = mom_in[momIdx<MomentId::rho>(index)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = ((4.f + 3.f * myyI) * rhoI) / 3.f * (1.f + uy);
    mom_out[momIdx<MomentId::mxx>(index)] = (18.f * mxxI * (1.f + uy)) / (20.f + 15.f * myyI);
    mom_out[momIdx<MomentId::mxy>(index)] = (-18.f * mxyI * (1.f + uy) + 4.f * ux + 3.f * myyI * ux) / (-12.f - 9.f * myyI);
    mom_out[momIdx<MomentId::myy>(index)] = (2.f + 15.f * myyI - 6.f * uy + 9.f * myyI * uy) / (12.f + 9.f * myyI);
}

__device__ void inlet_south(CInt *I_s, int x, int y,
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

    float rho = mom_in[momIdx<MomentId::rho>(index)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = ((-4.f - 3.f * myyI) * rhoI) / 3.f * (-1.f + uy);
    mom_out[momIdx<MomentId::mxx>(index)] = (18.f * mxxI * (-1.f + uy)) / (20.f + 15.f * myyI);
    mom_out[momIdx<MomentId::mxy>(index)] = (18.f * mxyI * (1.f - uy) + 4.f * ux + 3.f * myyI * ux) / (12.f + 9.f * myyI);
    mom_out[momIdx<MomentId::myy>(index)] = (2.f + 15.f * myyI + 6.f * uy - 9.f * myyI * uy) / (12.f + 9.f * myyI);
}

__device__ void inlet_east(CInt *I_s, int x, int y,
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

    float rho = mom_in[momIdx<MomentId::rho>(index)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = ((4.f + 3.f * mxxI) * rhoI) / 3.f * (1.f + ux);
    mom_out[momIdx<MomentId::mxx>(index)] = (2.f + 15.f * mxxI - 6.f * ux + 9.f * mxxI * ux) / (12.f + 9.f * mxxI);
    mom_out[momIdx<MomentId::mxy>(index)] = (18.f * mxyI * (1.f + ux) - 4.f * uy - 3.f * mxxI * uy) / (12.f + 9.f * mxxI);
    mom_out[momIdx<MomentId::myy>(index)] = (18.f * myyI * (1.f + ux)) / (20.f + 15.f * myyI);
}

__device__ void inlet_west(CInt *I_s, int x, int y,
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

    float rho = mom_in[momIdx<MomentId::rho>(index)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::rho>(index)] = ((-4.f - 3.f * mxxI) * rhoI) / 3.f * (-1.f + ux);
    mom_out[momIdx<MomentId::mxx>(index)] = (2.f + 15.f * mxxI + 6.f * ux - 9.f * mxxI * ux) / (12.f + 9.f * mxxI);
    mom_out[momIdx<MomentId::mxy>(index)] = (18.f * mxyI * (1.f - ux) + 4.f * uy + 3.f * mxxI * uy) / (12.f + 9.f * mxxI);
    mom_out[momIdx<MomentId::myy>(index)] = (18.f * myyI * (-1.f - ux)) / (20.f + 15.f * myyI);
}