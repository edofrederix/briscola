#ifndef DIRK2_H
#define DIRK2_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class DIRK2
:
    public RungeKuttaScheme
{
    // Butcher tableau
    List<scalarList> a_;

public:

    // Runtime type information
    TypeName("DIRK2");

    // Construct from mesh
    DIRK2(const fvMesh& fvMsh);

    // Destructor
    ~DIRK2()
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
    #include "DIRK2.C"
#endif

#endif
