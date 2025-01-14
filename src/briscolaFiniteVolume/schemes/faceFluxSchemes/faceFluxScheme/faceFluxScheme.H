#ifndef faceFluxScheme_H
#define faceFluxScheme_H

#include "scheme.H"
#include "FaceSpace.H"
#include "meshFields.H"
#include "meshLevels.H"
#include "meshDirections.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class faceFluxScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("faceFluxScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        faceFluxScheme,
        Istream,
        (
            const fvMesh& fvMsh,
            Istream& is
        ),
        (fvMsh, is)
    );

    // Construct from stream and mesh
    faceFluxScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Copy constructor
    faceFluxScheme(const faceFluxScheme&);

    // Selector
    static autoPtr<faceFluxScheme> New(const fvMesh&, const word);

    // Destructor
    virtual ~faceFluxScheme();

    // Face flux of a colocated vector field
    virtual tmp<colocatedFaceScalarField> faceFlux
    (
        const colocatedVectorField&
    ) = 0;

    // Face flux of a colocated tensor field
    virtual tmp<colocatedFaceVectorField> faceFlux
    (
        const colocatedTensorField&
    ) = 0;

    // Face flux of a staggered scalar field
    virtual tmp<staggeredFaceScalarField> faceFlux
    (
        const staggeredScalarField&
    ) = 0;
};

// Reconstruct colocated flux field from staggered velocity

tmp<colocatedFaceScalarField> coloFaceFlux
(
    const staggeredScalarField&
);

tmp<colocatedFaceScalarField> coloFaceFlux
(
    const tmp<staggeredScalarField>&
);

}

}

}

#endif
