#ifndef Ascher121_H
#define Ascher121_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher121
:
    public RungeKuttaScheme
{
    // Butcher tableaus
    scalarSquareMatrix a_;
    scalarSquareMatrix b_;

public:

    // Runtime type information
    TypeName("Ascher121");

    // Construct from mesh
    Ascher121(const fvMesh& fvMsh);

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
    #include "Ascher121.C"
#endif

#endif
