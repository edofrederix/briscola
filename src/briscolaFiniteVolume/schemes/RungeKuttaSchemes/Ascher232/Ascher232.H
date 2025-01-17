#ifndef Ascher232_H
#define Ascher232_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher232
:
    public RungeKuttaScheme
{
    // Butcher tableaus
    List<scalarList> a_;
    List<scalarList> b_;

public:

    // Runtime type information
    TypeName("Ascher232");

    // Construct from mesh
    Ascher232(const fvMesh& fvMsh);

    // Destructor
    ~Ascher232()
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
    #include "Ascher232.C"
#endif

#endif
