#ifndef vof_H
#define vof_H

#include "fvMesh.H"
#include "meshFields.H"
#include "regIOobject.H"

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
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the IO dictionary
    const IOdictionary& dict_;

    // Volume fraction field
    colocatedScalarField alpha_;

private:

    // Disallow default bitwise assignment
    void operator=(const vof&);

public:

    // Void fraction threshold
    static const scalar threshold;

    // Runtime type information
    TypeName("vof");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        vof,
        dictionary,
        (const IOdictionary& dict, const fvMesh& fvMsh),
        (dict, fvMsh)
    );

    // Construct from IOdictionary
    vof(const IOdictionary&, const fvMesh&);

    // Copy constructor
    vof(const vof&);

    // Clone vof
    virtual autoPtr<vof> clone() const = 0;

    // Selector
    static autoPtr<vof> New(const IOdictionary&, const fvMesh&);

    // Destructor
    virtual ~vof();

    // Needed by regIOobject
    virtual bool writeData(Ostream&) const { return true; }

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return access to the volume fraction field
    inline colocatedScalarField& alpha() { return alpha_; }

    // Return const access to the volume fraction field
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Solve the Vof advection equation given a flux field
    virtual void solve(const colocatedFaceScalarField&) = 0;
};

}

}

}

#endif
