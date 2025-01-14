#ifndef Ascher121_H
#define Ascher121_H

#include "RungeKuttaScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace RungeKuttaSchemes
{

class Ascher121
:
    public RungeKuttaScheme
{
    // Butcher tableaus
    List<scalarList> a_;
    List<scalarList> b_;

public:

    // Runtime type information
    TypeName("Ascher121");

    // Construct from mesh
    Ascher121(const fvMesh& fvMsh);

    // Destructor
    ~Ascher121()
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
    #include "Ascher121.C"
#endif

#endif
