#ifndef surfaceTensionScheme_H
#define surfaceTensionScheme_H

#include "dictionary.H"
#include "colocatedFields.H"
#include "normalScheme.H"
#include "curvatureScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class twoPhaseModel;

class surfaceTensionScheme
:
    public colocatedVectorField
{
protected:

    // Reference to the two-phase model
    const twoPhaseModel& tpm_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the dictionary
    const dictionary& dict_;

    // Reference to the normal scheme
    const normalScheme& normal_;

    // Reference to the void fraction
    const colocatedScalarField& alpha_;

    // Staggered surface tension force field, only created on structured meshes
    autoPtr<staggeredScalarField> stagForcePtr_;

    // Surface tension field
    colocatedScalarField sigma_;

    // Temporal fix
    scalar sigmaValue_;

    // Curvature scheme
    autoPtr<curvatureScheme> curvatureSchemePtr_;

    // Surface Tension Potential
    colocatedLowerFaceScalarField surfaceTensionPotential_;

public:

    // Runtime type information
    TypeName("surfaceTensionScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        surfaceTensionScheme,
        dictionary,
        (
            const twoPhaseModel& tpm,
            const dictionary& dict
        ),
        (tpm, dict)
    );

    // Construct from two-phase model
    surfaceTensionScheme
    (
        const twoPhaseModel&,
        const dictionary&
    );

    // Copy constructor
    surfaceTensionScheme(const surfaceTensionScheme&);

    // Clone surfaceTensionScheme
    virtual autoPtr<surfaceTensionScheme> clone() const = 0;

    // Selector from two-phase model
    static autoPtr<surfaceTensionScheme> New
    (
        const twoPhaseModel&,
        const dictionary&
    );

    // Destructor
    virtual ~surfaceTensionScheme();

    // Access to the two-phase model
    inline const twoPhaseModel& tpm() const { return tpm_; }

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the normal scheme
    inline const normalScheme& normal() const { return normal_; }

    // Access to the volume fraction field
    inline const colocatedScalarField& alpha() const { return alpha_; }

    // Access to the colocated surface tension force
    inline const colocatedVectorField& coloForce() const
    {
        return *this;
    }

    // Access to the staggered surface tension force
    inline const staggeredScalarField& stagForce() const
    {
        return stagForcePtr_();
    }

    // Access to the curvature scheme
    inline const curvatureScheme& kappa() const
    {
        return curvatureSchemePtr_();
    }

    // Access to the sigma field
    inline const colocatedScalarField& sigma() const { return sigma_; }

    // Access sigma value
    inline scalar sigmaValue() {return sigmaValue_; }

    // Correct the surface tension model
    virtual void correct();

    // Access to the surface tension porential
    inline const colocatedLowerFaceScalarField& potential() const
    {
        return surfaceTensionPotential_;
    }
};

}

}

}

#endif
