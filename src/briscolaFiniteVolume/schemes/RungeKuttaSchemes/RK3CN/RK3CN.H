#ifndef RK3CN_H
#define RK3CN_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class RK3CN
:
    public RK3
{
    // Butcher tableau
    List<scalarList> b_;

public:

    // Runtime type information
    TypeName("RK3CN");

    // Construct from mesh
    RK3CN(const fvMesh& fvMsh);

    // Destructor
    ~RK3CN()
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
    #include "RK3CN.C"
#endif

#endif
