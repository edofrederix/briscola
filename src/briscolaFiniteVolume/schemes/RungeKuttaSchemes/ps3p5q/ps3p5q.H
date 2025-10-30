#ifndef ps3p5q_H
#define ps3p5q_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class ps3p5q
:
    public RungeKuttaScheme
{
    // Butcher tableau
    scalarSquareMatrix a_;

public:

    // Runtime type information
    TypeName("ps3p5q");

    // Construct from mesh
    ps3p5q(const fvMesh& fvMsh);

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
    #include "ps3p5q.C"
#endif

#endif
