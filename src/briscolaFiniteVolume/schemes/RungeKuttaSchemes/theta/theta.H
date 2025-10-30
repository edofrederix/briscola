#ifndef theta_H
#define theta_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class theta
:
    public RungeKuttaScheme
{
    // Butcher tableau
    scalarSquareMatrix a_;

public:

    // Runtime type information
    TypeName("theta");

    // Construct from mesh
    theta(const fvMesh& fvMsh);

    // Construct from mesh and theta value
    theta(const fvMesh& fvMsh, const scalar theta);

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
    #include "theta.C"
#endif

#endif
