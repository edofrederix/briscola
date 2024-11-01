#ifndef Ascher233_H
#define Ascher233_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher233
:
    public RungeKuttaScheme
{
    // Butcher tableaus
    List<scalarList> a_;
    List<scalarList> b_;

public:

    // Runtime type information
    TypeName("Ascher233");

    // Construct from mesh
    Ascher233(const fvMesh& fvMsh);

    // Destructor
    ~Ascher233()
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
    #include "Ascher233.C"
#endif

#endif
