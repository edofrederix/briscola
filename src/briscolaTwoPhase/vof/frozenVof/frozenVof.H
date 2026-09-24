#ifndef frozenVof_H
#define frozenVof_H

#include "vof.H"
#include "LVE.H"
#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class frozenVof
:
    public vof
{
public:

    // Runtime type information
    TypeName("frozen");

    // Construct from components
    frozenVof
    (
        const fvMesh&,
        const dictionary&,
        normalScheme&,
        vofField&
    );

    // Copy constructor
    frozenVof(const frozenVof&);

    // Clone frozenVof
    virtual autoPtr<vof> clone() const
    {
        return autoPtr<vof>(new frozenVof(*this));
    }

    // Solve the Vof advection equation given a flux
    void solve(const colocatedScalarFaceField&);
};

}

}

}

#endif
