#pragma once

struct moments
{
    float *mom;
    float *mom_host;

    float *layer[Geometry::LNY];
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