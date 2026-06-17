#include "propagation.cuh"

__device__ void propagation(int x, int y, float **mom, float *layer, uint8_t *mask, uint8_t *node, int step_i)
{
   applyBoundary(x, y, node, mask, mom[step_i], mom[1 - step_i]);
}

__device__ void applyBoundary(int x, int y, uint8_t *node, uint8_t *mask, float *mom_in,
                              float *mom_out)
{
   const int index = grid_id();

   if (node[index] & to_u8(domainTags::Boundary))
   {
      boundary(x, y, mask[index], node[index], mom_in,
               mom_out);
   }
   else
   {
      center(x, y, mom_in,
             mom_out);
   }
}

__device__ void center(int x, int y,
                       float *mom_in,
                       float *mom_out)
{
   int index = grid_id();
   float rho = 0.f;
   float ux = 0.f;
   float uy = 0.f;
   float mxx = 0.f;
   float mxy = 0.f;
   float myy = 0.f;

#pragma unroll
   for (int i = 0; i < Q; i++)
   {
      int index_from = from_id(x, y, i);

      float fi = f_i(index_from, i, mom_in);

      rho += fi;

      ux += fi * (float)c_ix[i];
      uy += fi * (float)c_iy[i];

      mxx += fi * (c_ix[i] * c_ix[i] - inv_as2);
      mxy += fi * (c_ix[i] * c_iy[i]);
      myy += fi * (c_iy[i] * c_iy[i] - inv_as2);
   }
   mom_out[momIdx<MomentId::rho>(index)] = rho;

   mom_out[momIdx<MomentId::ux>(index)] = ux / rho;
   mom_out[momIdx<MomentId::uy>(index)] = uy / rho;

   mom_out[momIdx<MomentId::mxx>(index)] = mxx / rho;
   mom_out[momIdx<MomentId::myy>(index)] = myy / rho;
   mom_out[momIdx<MomentId::mxy>(index)] = mxy / rho;
}

__device__ void boundary(int x, int y, uint8_t mask_in, uint8_t node_in, float *mom_in,
                         float *mom_out)
{

   int index = grid_id();
   float sum_fi = 0.f;
   float mxy_I = 0.f;

   for (int i = 0; i < Q; i++)
   {
      if (mask_in & (1u << i))
      {
         int index_from = from_id(x, y, i);

         float fi = f_i(index_from, i, mom_in);
         sum_fi += fi;
         mxy_I += fi * c_ix[i] * c_iy[i];
      }
   }

   mxy_I /= sum_fi;

   switch (static_cast<Boundary>(node_in))
   {
   case Boundary::East:
      Constants<Boundary::East> c(mask_in);
      mom_out[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
      mom_out[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * mom_out[momIdx<MomentId::mxy>(index)]);
      break;
   case Boundary::North:
      Constants<Boundary::North> c(mask_in);
      mom_out[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
      mom_out[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * mom_out[momIdx<MomentId::mxy>(index)]);
      break;
   case Boundary::West:
      Constants<Boundary::West> c(mask_in);
      mom_out[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
      mom_out[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * mom_out[momIdx<MomentId::mxy>(index)]);
      break;
   case Boundary::South:
      Constants<Boundary::South> c(mask_in);
      mom_out[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
      mom_out[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * mom_out[momIdx<MomentId::mxy>(index)]);
      break;
   case Boundary::Northeast:
      Constants<Boundary::Northeast> c(mask_in);
      mom_out[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
      mom_out[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * mom_out[momIdx<MomentId::mxy>(index)]);
      break;
   case Boundary::Northwest:
      Constants<Boundary::Northwest> c(mask_in);
      mom_out[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
      mom_out[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * mom_out[momIdx<MomentId::mxy>(index)]);
      break;
   case Boundary::Southwest:
      Constants<Boundary::Southwest> c(mask_in);
      mom_out[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
      mom_out[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * mom_out[momIdx<MomentId::mxy>(index)]);
      break;
   case Boundary::Southeast:
      Constants<Boundary::Southeast> c(mask_in);
      mom_out[momIdx<MomentId::mxy>(index)] = (c.C3 - c.C1 * mxy_I) / (c.C2 * mxy_I - c.C4);
      mom_out[momIdx<MomentId::rho>(index)] = sum_fi / (c.C1 + c.C2 * mom_out[momIdx<MomentId::mxy>(index)]);
      break;
   default:
      break;
   }
}
