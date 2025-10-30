#ifndef Ascher122_H
#define Ascher122_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher122
:
    public RungeKuttaScheme
{
    // Butcher tableaus
    scalarSquareMatrix a_;
    scalarSquareMatrix b_;

public:

    // Runtime type information
    TypeName("Ascher122");

    // Construct from mesh
    Ascher122(const fvMesh& fvMsh);

    // Return the Butcher tableau for source A
    virtual const scalarSquareMatrix& a() const
    {
        return a_;
    }

    // Return the Butcher tableau for source B
    virtual const scalarSquareMatrix& b() const
    {
        return b_;
    }
};

}

}

}

}

#ifdef NoRepository
    #include "Ascher122.C"
#endif

#endif
