#ifndef Ascher111_H
#define Ascher111_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher111
:
    public RungeKuttaScheme
{
    // Butcher tableaus
    List<scalarList> a_;
    List<scalarList> b_;

public:

    // Runtime type information
    TypeName("Ascher111");

    // Construct from mesh
    Ascher111(const fvMesh& fvMsh);

    // Destructor
    ~Ascher111()
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
    #include "Ascher111.C"
#endif

#endif
