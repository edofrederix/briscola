#ifndef midPoint_H
#define midPoint_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class midPoint
:
    public RungeKuttaScheme
{
    // Butcher tableau
    scalarSquareMatrix a_;

public:

    // Runtime type information
    TypeName("midPoint");

    // Construct from mesh
    midPoint(const fvMesh& fvMsh);

    // Return the Butcher tableau for source A
    virtual const scalarSquareMatrix& a() const
    {
        return a_;
    }

    // Return the Butcher tableau for source B
    virtual const scalarSquareMatrix& b() const
    {
        return a_;
    }
};

}

}

}

}

#ifdef NoRepository
    #include "midPoint.C"
#endif

#endif
