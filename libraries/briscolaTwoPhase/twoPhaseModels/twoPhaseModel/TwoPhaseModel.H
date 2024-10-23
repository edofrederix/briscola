#ifndef TwoPhaseModel_H
#define TwoPhaseModel_H

#include "twoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class MeshType>
class TwoPhaseModel
:
    public twoPhaseModel
{
protected:

    // Mixture mass density
    meshField<scalar,MeshType> rho_;

    // Mixture dynamic viscosity
    meshField<lowerFaceScalar,MeshType> mu_;

private:

    // Disallow default bitwise assignment
    void operator=(const TwoPhaseModel&);

public:

    // Type of mesh
    typedef MeshType meshType;

    // Construct from IOdictionary
    TwoPhaseModel(const fvMesh&, const IOdictionary&);

    // Copy constructor
    TwoPhaseModel(const TwoPhaseModel&);

    // Destructor
    virtual ~TwoPhaseModel();

    // Selector returning templated type
    static autoPtr<TwoPhaseModel<MeshType>> New
    (
        const fvMesh&,
        const IOdictionary&
    );

    // Return access to the mixture mass density
    meshField<scalar,MeshType>& rho() { return rho_; }

    // Return const access to the mixture mass density
    const meshField<scalar,MeshType>& rho() const { return rho_; }

    // Return access to the mixture dynamic viscosity
    meshField<lowerFaceScalar,MeshType>& mu() { return mu_; }

    // Return const access to the mixture dynamic viscosity
    const meshField<lowerFaceScalar,MeshType>& mu() const { return mu_; }

    // Return an interpolated value of the volume fraction
    virtual tmp<meshField<lowerFaceScalar,MeshType>> faceAlpha() const;

    // Return the volume fraction for the templated mesh type
    virtual tmp<meshField<scalar,MeshType>> meshAlpha() const;

    // Return the colocated face flux
    virtual tmp<colocatedLowerFaceScalarField> coloFaceFlux() const;

    // Return the gravitational acceleration
    virtual List<typename MeshType::vectorType> gravity() const;
};

}

}

}

#ifdef NoRepository
    #include "TwoPhaseModel.C"
#endif

#endif
