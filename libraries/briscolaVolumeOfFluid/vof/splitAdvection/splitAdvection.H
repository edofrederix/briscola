#ifndef splitAdvection_H
#define splitAdvection_H

#include "vof.H"
#include "LVE.H"
#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class splitAdvection
:
    public vof
{
protected:

    // Flux field
    colocatedFaceScalarField flux_;

    // Rectilinear bool (stored here for efficiency)
    bool rectilinear_;

    // Local volume enforcement problem solver
    LVE lve_;

    // Normal scheme
    autoPtr<normalScheme> normalSchemePtr_;

    // Reset the flux
    inline void resetFlux()
    {
        flux_ = Zero;
    }

    // Update the flux in given direction
    void updateFlux
    (
        const colocatedFaceScalarField&,
        const label
    );

private:

    // Disallow default bitwise assignment
    void operator=(const splitAdvection&);

public:

    // Runtime type information
    TypeName("splitAdvection");

    // Construct from IOdictionary
    splitAdvection(const IOdictionary&, const fvMesh&);

    // Copy constructor
    splitAdvection(const splitAdvection&);

    // Clone splitAdvection
    virtual autoPtr<vof> clone() const
    {
        return autoPtr<vof>(new splitAdvection(*this));
    }

    // Destructor
    ~splitAdvection();

    // Reference to the LVE solver
    inline const LVE& lve() const { return lve_; }

    // Reference to the normal scheme
    inline normalScheme& normal() { return normalSchemePtr_(); }

    // Solve the Vof advection equation given a flux
    void solve(const colocatedFaceScalarField&);
};

}

}

}

#endif
