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

    float rho = mom_in[momIdx<MomentId::rho>(index)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::mxx>(index)] = (6.f * mxxI * rhoI) / (5.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = (-6.f * mxyI * rhoI + rho * ux) / (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (-rho - 9.f * myyI * rhoI - 3.f * rho * uy) / (6.f * rho);
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

    float rho = mom_in[momIdx<MomentId::rho>(index)];
    float mxx = mom_in[momIdx<MomentId::mxx>(index)];
    float mxy = mom_in[momIdx<MomentId::mxy>(index)];
    float myy = mom_in[momIdx<MomentId::myy>(index)];

    mom_out[momIdx<MomentId::mxx>(index)] = (6.f * mxxI * rhoI) / (5.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = (-6.f * mxyI * rhoI - rho * ux) / (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (-rho - 9.f * myyI * rhoI - 3.f * rho * uy) / (6.f * rho);
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

    mom_out[momIdx<MomentId::mxx>(index)] = (-rho - 9.f * mxxI * rhoI + 3.f * rho * ux) / (6.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = (-6.f * mxyI * rhoI + rho * uy) / (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (6.f * myyI * rhoI) / (5.f * rho);
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

    mom_out[momIdx<MomentId::mxx>(index)] = (-rho - 9.f * mxxI * rhoI - 3.f * rho * ux) / (6.f * rho);
    mom_out[momIdx<MomentId::mxy>(index)] = (-6.f * mxyI * rhoI - rho * uy) / (3.f * rho);
    mom_out[momIdx<MomentId::myy>(index)] = (6.f * myyI * rhoI) / (5.f * rho);
}