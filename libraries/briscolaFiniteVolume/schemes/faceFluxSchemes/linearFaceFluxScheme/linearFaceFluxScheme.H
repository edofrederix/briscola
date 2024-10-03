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
            LowerFaceSpace<typename innerProduct<Type,vector>::type>,
            colocated
        >
    > coloFaceFlux
    (
        const meshField<Type,colocated>&
    );

public:

    // Runtime type information
    TypeName("linear");

    // Construct from dictionary and mesh
    linearFaceFluxScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    linearFaceFluxScheme(const fvMesh&);

    // Face flux of a colocated scalar field
    virtual tmp<colocatedLowerFaceScalarField> faceFlux
    (
        const colocatedVectorField&
    );

    // Face flux of a colocated tensor field
    virtual tmp<colocatedLowerFaceVectorField> faceFlux
    (
        const colocatedTensorField&
    );

    // Face flux of a staggered scalar field
    virtual tmp<staggeredLowerFaceScalarField> faceFlux
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
