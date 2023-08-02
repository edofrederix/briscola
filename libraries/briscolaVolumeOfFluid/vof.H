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

    // Local volume enforcement problem solver
    LVE lve_;

    // Normal scheme
    autoPtr<normalScheme> normalSchemePtr_;

    // Disallow default bitwise assignment
    void operator=(const vof&);

    // Disallow copy constructor
    vof(const vof&);

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
    inline const normalScheme& normal() const { return normalSchemePtr_(); }
};

}

}

}

#endif
