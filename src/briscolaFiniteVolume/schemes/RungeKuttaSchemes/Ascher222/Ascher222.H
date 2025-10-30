#ifndef Ascher222_H
#define Ascher222_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher222
:
    public RungeKuttaScheme
{
    // Butcher tableaus
    scalarSquareMatrix a_;
    scalarSquareMatrix b_;

public:

    // Runtime type information
    TypeName("Ascher222");

    // Construct from mesh
    Ascher222(const fvMesh& fvMsh);

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
    #include "Ascher222.C"
#endif

#endif
