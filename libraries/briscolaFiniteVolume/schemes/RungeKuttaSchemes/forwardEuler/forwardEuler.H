#ifndef forwardEuler_H
#define forwardEuler_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class forwardEuler
:
    public RungeKuttaScheme
{
    // Butcher tableau
    List<scalarList> a_;

public:

    // Runtime type information
    TypeName("forwardEuler");

    // Construct from mesh
    forwardEuler(const fvMesh& fvMsh);

    // Destructor
    ~forwardEuler()
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
    #include "forwardEuler.C"
#endif

#endif
