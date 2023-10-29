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
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Reference to the dictionary
    const dictionary& dict_;

    // Reference to the normal scheme
    const normalScheme& normal_;

    // Reference to the void fraction
    const colocatedScalarField& alpha_;

    // Colocated surface tension force field
    colocatedVectorField coloForce_;

    // Staggered surface tension force field, only created on structured meshes
    autoPtr<staggeredScalarField> stagForcePtr_;

    // Surface tension field
    colocatedScalarField sigma_;

    // Curvature scheme
    autoPtr<curvatureScheme> curvatureSchemePtr_;

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
            const fvMesh& fvMsh,
            const dictionary& dict,
            const normalScheme& normal,
            const colocatedScalarField& alpha
        ),
        (fvMsh, dict, normal, alpha)
    );

    // Construct from components
    surfaceTensionScheme
    (
        const fvMesh&,
        const dictionary&,
        const normalScheme&,
        const colocatedScalarField&
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

    // Selector from components
    static autoPtr<surfaceTensionScheme> New
    (
        const fvMesh&,
        const dictionary&,
        const normalScheme&,
        const colocatedScalarField&
    );

    // Destructor
    virtual ~surfaceTensionScheme();

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
        return coloForce_;
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

    // Correct the surface tension model
    virtual void correct() = 0;
};

}

}

}

#endif
