#ifndef Alexander3_H
#define Alexander3_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Alexander3
:
    public RungeKuttaScheme
{
    // Butcher tableau
    scalarSquareMatrix a_;

public:

    // Runtime type information
    TypeName("Alexander3");

    // Construct from mesh
    Alexander3(const fvMesh& fvMsh);

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
    #include "Alexander3.C"
#endif

#endif
