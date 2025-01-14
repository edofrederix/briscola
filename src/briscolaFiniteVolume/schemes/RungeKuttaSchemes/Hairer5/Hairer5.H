#ifndef Hairer5_H
#define Hairer5_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Hairer5
:
    public RungeKuttaScheme
{
    // Butcher tableau
    List<scalarList> a_;

public:

    // Runtime type information
    TypeName("Hairer5");

    // Construct from mesh
    Hairer5(const fvMesh& fvMsh);

    // Destructor
    ~Hairer5()
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
    #include "Hairer5.C"
#endif

#endif
