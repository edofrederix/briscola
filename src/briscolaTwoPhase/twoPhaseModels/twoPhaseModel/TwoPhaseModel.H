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
    faceField<scalar,MeshType> mu_;

    // Cached g*h*A face field
    autoPtr<colocatedScalarFaceField> ghAPtr_;

    // Function to set the g*h*A field
    void setGhA()
    {
        if (!ghAPtr_.valid())
        {
            const colocatedVectorFaceField& fc =
                fvMsh_.template metrics<colocated>().faceCenters();

            const colocatedScalarFaceField& fa =
                fvMsh_.template metrics<colocated>().faceAreas();

            tmp<colocatedScalarFaceField> tghA(fa);
            colocatedScalarFaceField& ghA = tghA.ref();

            forAllFaces(ghA, fd, i, j, k)
                ghA[fd](i,j,k) *= -(this->g() & fc[fd](i,j,k));

            ghAPtr_.set(tghA.ptr());
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
    faceField<scalar,MeshType>& mu() { return mu_; }

    // Return const access to the mixture dynamic viscosity
    const faceField<scalar,MeshType>& mu() const { return mu_; }

    // Return an interpolated value of the volume fraction
    virtual tmp<faceField<scalar,MeshType>> faceAlpha() const;

    // Return the volume fraction for the templated mesh type
    virtual tmp<meshField<scalar,MeshType>> meshAlpha() const;

    // Return the colocated face flux
    virtual tmp<colocatedScalarFaceField> coloFaceFlux() const;

    // Return the colocated mass density
    virtual tmp<colocatedScalarField> coloRho() const = 0;

    // Compute the mean mixture mass density
    virtual scalar rhoMean() const;

    // Return the buoyancy term defined as (rho - rhoMean)*g
    virtual tmp<meshField<MeshVectorType,MeshType>> buoyancy() const;

    // Return the interfacial flux
    tmp<colocatedScalarFaceField> flux();
};

}

}

}

#endif
