#ifndef TwoPhaseModel_H
#define TwoPhaseModel_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Mesh-dependent part of the two-phase model

template<class MeshType>
class TwoPhaseModel
:
    public twoPhaseModel
{
protected:

    // Face volume fraction
    faceField<scalar,MeshType> faceAlpha_;

    // Mixture dynamic viscosity
    faceField<scalar,MeshType> mu_;

    // Return the mesh-dependent colocated face flux
    virtual tmp<colocatedScalarFaceField> coloFaceFlux() const;

    // Correct the face volume fraction
    virtual void correctFaceAlpha();

public:

    // Type of mesh
    typedef MeshType meshType;

    // Type of vector on this mesh type
    typedef typename MeshType::vectorType MeshVectorType;

    // Construct from IOdictionary
    TwoPhaseModel(const fvMesh&, const IOdictionary&);

    // Copy constructor
    TwoPhaseModel(const TwoPhaseModel&);

    // Selector returning templated type
    static autoPtr<TwoPhaseModel<MeshType>> New
    (
        const fvMesh&,
        const IOdictionary&
    );

    // Return access to the face volume fraction
    faceField<scalar,MeshType>& faceAlpha() { return faceAlpha_; }

    // Return const access to the face volume fraction
    const faceField<scalar,MeshType>& faceAlpha() const { return faceAlpha_; }

    // Return access to the mixture dynamic viscosity
    faceField<scalar,MeshType>& mu() { return mu_; }

    // Return const access to the mixture dynamic viscosity
    const faceField<scalar,MeshType>& mu() const { return mu_; }

    // Return the colocated mass density
    virtual tmp<colocatedScalarField> coloRho() const;

    // Return the staggered mass density (not defined for colocated mesh type)
    virtual tmp<staggeredScalarField> stagRho() const;

    // Return the buoyancy term defined as (rho - rhoMean)*g
    virtual tmp<meshField<MeshVectorType,MeshType>> buoyancy() const;
};

}

}

}

#endif
