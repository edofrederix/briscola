#ifndef Ascher122_H
#define Ascher122_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher122
:
    public RungeKuttaScheme
{
    // Butcher tableaus
    List<scalarList> a_;
    List<scalarList> b_;

public:

    // Runtime type information
    TypeName("Ascher122");

    // Construct from mesh
    Ascher122(const fvMesh& fvMsh);

    // Destructor
    ~Ascher122()
    {}

    // Return the Butcher tableau for source A
    virtual const List<scalarList>& a() const
    {
        return a_;
    }

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
    #include "Ascher122.C"
#endif

#endif
