#ifndef Ascher232_H
#define Ascher232_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher232
:
    public RungeKuttaScheme
{
    // Butcher tableaus
    scalarSquareMatrix a_;
    scalarSquareMatrix b_;

public:

    // Runtime type information
    TypeName("Ascher232");

    // Construct from mesh
    Ascher232(const fvMesh& fvMsh);

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
    #include "Ascher232.C"
#endif

#endif
