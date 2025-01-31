#ifndef RK4CN_H
#define RK4CN_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class RK4CN
:
    public RK4
{
    // Butcher tableau
    List<scalarList> b_;

public:

    // Runtime type information
    TypeName("RK4CN");

    // Construct from mesh
    RK4CN(const fvMesh& fvMsh);

    // Destructor
    ~RK4CN()
    {}

    // Return the Butcher tableau for source B
    virtual const List<scalarList>& b() const
    {
        return b_;
    }
};

}

}

}

}

#ifdef NoRepository
    #include "RK4CN.C"
#endif

#endif
