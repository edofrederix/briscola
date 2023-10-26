#ifndef twoPhaseModel_H
#define twoPhaseModel_H

#include "fvMesh.H"
#include "meshFields.H"
#include "regIOobject.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class twoPhaseModel
:
    public regIOobject
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the IO dictionary
    const IOdictionary& dict_;

    // Colocated mixture mass density and dynamic viscosity
    colocatedScalarField rhoc_;
    colocatedFaceScalarField muc_;

    // Staggered mixture mass density and dynamic viscosity, only created on
    // structured meshes
    autoPtr<staggeredScalarField> rhosPtr_;
    autoPtr<staggeredFaceScalarField> musPtr_;

    // Function to set the density and viscosity restriction schemes. Can be
    // overloaded by derived types.
    virtual void setRestrictionSchemes();

private:

    // Disallow default bitwise assignment
    void operator=(const twoPhaseModel&);

public:

    // Void fraction threshold
    static const scalar threshold;

    // Runtime type information
    TypeName("twoPhaseModel");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        twoPhaseModel,
        dictionary,
        (const IOdictionary& dict, const fvMesh& fvMsh),
        (dict, fvMsh)
    );

    // Construct from IOdictionary
    twoPhaseModel(const IOdictionary&, const fvMesh&);

    // Copy constructor
    twoPhaseModel(const twoPhaseModel&);

    // Clone twoPhaseModel
    virtual autoPtr<twoPhaseModel> clone() const = 0;

    // Selector
    static autoPtr<twoPhaseModel> New(const IOdictionary&, const fvMesh&);

    // Destructor
    virtual ~twoPhaseModel();

    // Needed by regIOobject
    virtual bool writeData(Ostream&) const { return true; }

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return const access to the colocated mixture mass density and dynamic
    // viscosity
    inline const colocatedScalarField& rhoc() const { return rhoc_; }
    inline const colocatedFaceScalarField& muc() const { return muc_; }

    // Return const access to the staggered mixture mass density and dynamic
    // viscosity fields (only valid on structured meshes)
    inline const staggeredScalarField& rhos() const { return rhosPtr_(); }
    inline const staggeredFaceScalarField& mus() const { return musPtr_(); }

    // Correct the two-phase model
    virtual void correct() = 0;
};

}

}

}

#endif
