#ifndef twoPhaseModel_H
#define twoPhaseModel_H

#include "fvMesh.H"
#include "meshFields.H"
#include "regIOobject.H"
#include "normalScheme.H"
#include "surfaceTensionScheme.H"

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

    // Volume fraction field
    colocatedScalarField alpha_;

    // Colocated phase densities
    colocatedScalarField rho1_;
    colocatedScalarField rho2_;

    // Colocated mixture mass density and dynamic viscosity
    colocatedScalarField rhoc_;
    colocatedFaceScalarField muc_;

    // Staggered phase densities
    autoPtr<staggeredScalarField> rho1Ptr_;
    autoPtr<staggeredScalarField> rho2Ptr_;

    // Staggered mixture mass density and dynamic viscosity, only created on
    // structured meshes
    autoPtr<staggeredScalarField> rhosPtr_;
    autoPtr<staggeredFaceScalarField> musPtr_;

    // Normal scheme
    autoPtr<normalScheme> normalSchemePtr_;

    // Surface tension scheme
    autoPtr<surfaceTensionScheme> surfaceTensionSchemePtr_;

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
        (
            const fvMesh& fvMsh,
            const IOdictionary& dict
        ),
        (fvMsh, dict)
    );

    // Construct from IOdictionary
    twoPhaseModel(const fvMesh&, const IOdictionary&);

    // Copy constructor
    twoPhaseModel(const twoPhaseModel&);

    // Clone twoPhaseModel
    virtual autoPtr<twoPhaseModel> clone() const = 0;

    // Selector
    static autoPtr<twoPhaseModel> New(const fvMesh&, const IOdictionary&);

    // Destructor
    virtual ~twoPhaseModel();

    // Check if cast to derived model is possible
    template<class Model>
    inline bool castable() const
    {
        return dynamic_cast<const Model*>(this);
    }

    // Cast to derived model, if possible. Returns error otherwise.
    template<class Model>
    inline Model& cast()
    {
        Model* tpmPtr = dynamic_cast<Model*>(this);

        if (!tpmPtr)
        {
            FatalErrorInFunction
                << "Two-phase model cannot be cast to " << Model::typeName
                << " because it is " << this->type()
                << abort(FatalError);
        }

        return *tpmPtr;
    }

    // Needed by regIOobject
    virtual bool writeData(Ostream&) const { return true; }

    // Return const access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return access to the volume fraction
    inline colocatedScalarField& alpha() { return alpha_; }

    // Return const access to the volume fraction
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Return const access to the mixture mass density
    template<class MeshType>
    const meshField<scalar,MeshType>& rho() const;

    // Return const access to the mixture dynamic viscosity
    template<class MeshType>
    const meshField<faceScalar,MeshType>& mu() const;

    // Return access the normal scheme
    normalScheme& normal() { return normalSchemePtr_(); }

    // Return const access to the normal scheme
    const normalScheme& normal() const { return normalSchemePtr_(); }

    // Return access the surface tension scheme
    surfaceTensionScheme& surfaceTension()
    {
        return surfaceTensionSchemePtr_();
    }

    // Return const access to the surface tension scheme
    const surfaceTensionScheme& surfaceTension() const
    {
        return surfaceTensionSchemePtr_();
    }

    // Return the mean mixture mass density defined as <rho> = (rho1+rho2)/2
    template<class MeshType>
    tmp<meshField<scalar,MeshType>> meanRho() const;

    // Correct the two-phase mixture
    virtual void correctMixture() = 0;

    // Correct the two-phase model
    virtual void correct() = 0;
};

}

}

}

#endif
