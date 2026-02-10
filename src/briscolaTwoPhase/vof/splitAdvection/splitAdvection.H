#ifndef splitAdvection_H
#define splitAdvection_H

#include "geometricVof.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class splitAdvection
:
    public geometricVof
{
protected:

    // Flux field
    colocatedScalarFaceField flux_;

    // Reset the flux
    inline void resetFlux()
    {
        flux_ = Zero;
    }

    // Update the flux in given direction
    void updateFlux
    (
        const colocatedScalarFaceField&,
        const label
    );

public:

    // Runtime type information
    TypeName("splitAdvection");

    // Construct from components
    splitAdvection
    (
        const fvMesh&,
        const dictionary&,
        normalScheme&,
        vofField&
    );

    // Copy constructor
    splitAdvection(const splitAdvection&);

    // Clone splitAdvection
    virtual autoPtr<vof> clone() const
    {
        return autoPtr<vof>(new splitAdvection(*this));
    }

    // Solve the Vof advection equation given a flux
    void solve(const colocatedScalarFaceField&);
};

}

}

}

#endif
