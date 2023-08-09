#ifndef vof_H
#define vof_H

#include "fvMesh.H"
#include "meshFields.H"
#include "regIOobject.H"
#include "Time.H"
#include "LVE.H"
#include "normalScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class vof
:
    public regIOobject
{
    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Volume fraction field
    colocatedScalarField alpha_;

    // Flux field
    colocatedFaceScalarField flux_;

    // Rectilinear bool (stored here for efficiency)
    bool rectilinear_;

    // Local volume enforcement problem solver
    LVE lve_;

    // Normal scheme
    autoPtr<normalScheme> normalSchemePtr_;

    // Void fraction threshold
    static const scalar threshold_;

    // Disallow default bitwise assignment
    void operator=(const vof&);

    // Disallow copy constructor
    vof(const vof&);

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

public:

    // Declare name of the class and its debug switch
    ClassName("vof");

    // Construct from IOdictionary
    vof(const IOdictionary&, const fvMesh&);

    // Destructor
    ~vof();

    // Needed by regIOobject
    virtual bool writeData(Ostream&) const { return true; }

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return access to the volume fraction field
    inline colocatedScalarField& alpha() { return alpha_; }

    // Return const access to the volume fraction field
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Reference to the LVE solver
    inline const LVE& lve() const { return lve_; }

    // Reference to the normal scheme
    inline normalScheme& normal() { return normalSchemePtr_(); }

    // Solve the Vof advection equation given a flux
    void solve(const colocatedFaceScalarField& phi);
};

}

}

}

#endif
