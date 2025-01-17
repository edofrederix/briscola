#ifndef linearFaceFluxScheme_H
#define linearFaceFluxScheme_H

#include "faceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class linearFaceFluxScheme
:
    public faceFluxScheme
{
    // Face flux of a colocated field
    template<class Type>
    tmp
    <
        meshField
        <
            FaceSpace<typename innerProduct<Type,vector>::type>,
            colocated
        >
    > coloFaceFlux
    (
        const meshField<Type,colocated>&
    );

public:

    // Runtime type information
    TypeName("linear");

    // Construct from stream and mesh
    linearFaceFluxScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Face flux of a colocated scalar field
    virtual tmp<colocatedFaceScalarField> faceFlux
    (
        const colocatedVectorField&
    );

    // Face flux of a colocated tensor field
    virtual tmp<colocatedFaceVectorField> faceFlux
    (
        const colocatedTensorField&
    );

    // Face flux of a staggered scalar field
    virtual tmp<staggeredFaceScalarField> faceFlux
    (
        const staggeredScalarField&
    );
};

}

}

}

#ifdef NoRepository
    #include "linearFaceFluxScheme.C"
#endif

#endif
