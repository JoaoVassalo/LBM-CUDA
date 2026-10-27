#pragma once

#include "simulationConfig.cuh"

struct moments
{
    varType *mom;
    varType *mom_host;
    varType *vtk;

    varType *layer[Geometry::LNY];
};

enum momId
{
    rho = 0,
    ux,
    uy,
    mxx,
    mxy,
    myy,
    count
};
