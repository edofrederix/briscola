#ifndef backwardEuler_H
#define backwardEuler_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class backwardEuler
:
    public RungeKuttaScheme
{
    // Butcher tableau
    List<scalarList> a_;

public:

    // Runtime type information
    TypeName("backwardEuler");

    // Construct from mesh
    backwardEuler(const fvMesh& fvMsh);

    // Destructor
    ~backwardEuler()
    {}

    // Return the Butcher tableau for source A
    virtual const List<scalarList>& a() const
    {
        return a_;
    }

    // Return the Butcher tableau for source B
    virtual const List<scalarList>& b() const
    {
        return a_;
    }
};

}

}

}

}

#ifdef NoRepository
    #include "backwardEuler.C"
#endif

#endif
