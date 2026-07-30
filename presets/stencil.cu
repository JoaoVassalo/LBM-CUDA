#include "stencil.cuh"

extern __constant__ varUnit w[Q] = {
    4.f / 9.f,
    1.f / 9.f, 1.f / 9.f, 1.f / 9.f, 1.f / 9.f,
    1.f / 36.f, 1.f / 36.f, 1.f / 36.f, 1.f / 36.f};
extern __constant__ varUnit c_ix[Q] = {0.f, 1.f, 0.f, -1.f, 0.f, 1.f, -1.f, -1.f, 1.f};
extern __constant__ varUnit c_iy[Q] = {0.f, 0.f, 1.f, 0.f, -1.f, 1.f, 1.f, -1.f, -1.f};

extern __constant__ int Is_SW[4] = {0, 3, 4, 7};
extern __constant__ int Os_SW[4] = {0, 1, 2, 5};

extern __constant__ int Is_SE[4] = {0, 1, 4, 8};
extern __constant__ int Os_SE[4] = {0, 2, 3, 6};

extern __constant__ int Is_NW[4] = {0, 2, 3, 6};
extern __constant__ int Os_NW[4] = {0, 1, 4, 8};

extern __constant__ int Is_NE[4] = {0, 1, 2, 5};
extern __constant__ int Os_NE[4] = {0, 3, 4, 7};

extern __constant__ int Is_N[6] = {0, 1, 2, 3, 5, 6};
extern __constant__ int Os_N[6] = {0, 1, 3, 4, 7, 8};

extern __constant__ int Is_S[6] = {0, 1, 3, 4, 7, 8};
extern __constant__ int Os_S[6] = {0, 1, 2, 3, 5, 6};

extern __constant__ int Is_W[6] = {0, 2, 3, 4, 6, 7};
extern __constant__ int Os_W[6] = {0, 1, 2, 4, 5, 8};

extern __constant__ int Is_E[6] = {0, 1, 2, 4, 5, 8};
extern __constant__ int Os_E[6] = {0, 2, 3, 4, 6, 7};