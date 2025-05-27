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
    meshField<faceScalar,MeshType> mu_;

    // Cached g*h*A face field
    autoPtr<colocatedFaceScalarField> ghAPtr_;

    // Function to set the g*h*A field
    void setGhA()
    {
        if (!ghAPtr_.valid())
        {
            const colocatedFaceVectorField& fc =
                fvMsh_.template metrics<colocated>().faceCenters();

            const colocatedFaceScalarField& fa =
                fvMsh_.template metrics<colocated>().faceAreas();

            ghAPtr_.set((- (this->g() & fc) * fa).ptr());
        }
    }

private:

    // Disallow default bitwise assignment
    void operator=(const TwoPhaseModel&);

public:

    // Type of mesh
    typedef MeshType meshType;

    // Type of vector on this mesh type
    typedef typename MeshType::vectorType MeshVectorType;

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
    meshField<faceScalar,MeshType>& mu() { return mu_; }

    // Return const access to the mixture dynamic viscosity
    const meshField<faceScalar,MeshType>& mu() const { return mu_; }

    // Return an interpolated value of the volume fraction
    virtual tmp<meshField<faceScalar,MeshType>> faceAlpha() const;

    // Return the volume fraction for the templated mesh type
    virtual tmp<meshField<scalar,MeshType>> meshAlpha() const;

    // Return the colocated face flux
    virtual tmp<colocatedFaceScalarField> coloFaceFlux() const;

    // Return the colocated mass density
    virtual tmp<colocatedScalarField> coloRho() const = 0;

    // Compute the mean mixture mass density
    virtual scalar rhoMean() const;

    // Return the buoyancy term defined as (rho - rhoMean)*g
    virtual tmp<meshField<MeshVectorType,MeshType>> buoyancy() const;

    // Return the interfacial flux
    tmp<colocatedFaceScalarField> flux();
};

}

}

}

#endif
