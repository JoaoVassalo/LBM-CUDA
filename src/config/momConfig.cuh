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

enum Boundary
{
    Center = 0,
    East,
    North,
    West,
    South,
    Northeast,
    Northwest,
    Southwest,
    Southeast
};

__host__ __device__ __forceinline__ Boundary defBC(int x, int y)
{
    if (x == 0 && y == 0)
    {
        return Boundary::Southwest;
    }
    else if (x == Geometry::NX - 1 && y == 0)
    {
        return Boundary::Southeast;
    }
    else if (x == 0 && y == Geometry::NY - 1)
    {
        return Boundary::Northwest;
    }
    else if (x == Geometry::NX - 1 && y == Geometry::NY - 1)
    {
        return Boundary::Northeast;
    }
    else if (y == 0)
    {
        return Boundary::South;
    }
    else if (y == Geometry::NY - 1)
    {
        return Boundary::North;
    }
    else if (x == 0)
    {
        return Boundary::West;
    }
    else if (x == Geometry::NX - 1)
    {
        return Boundary::East;
    }
    else
    {
        return Boundary::Center;
    }
}