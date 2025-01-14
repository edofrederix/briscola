#ifndef ps3p5q_H
#define ps3p5q_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class ps3p5q
:
    public RungeKuttaScheme
{
    // Butcher tableau
    List<scalarList> a_;

public:

    // Runtime type information
    TypeName("ps3p5q");

    // Construct from mesh
    ps3p5q(const fvMesh& fvMsh);

    // Destructor
    ~ps3p5q()
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
    #include "ps3p5q.C"
#endif

#endif
