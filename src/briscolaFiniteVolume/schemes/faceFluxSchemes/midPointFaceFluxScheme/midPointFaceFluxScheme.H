#ifndef midPointFaceFluxScheme_H
#define midPointFaceFluxScheme_H

#include "faceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class midPointFaceFluxScheme
:
    public faceFluxScheme
{
    // Face flux of a colocated field
    template<class Type>
    tmp
    <
        faceField
        <
            typename innerProduct<Type,vector>::type,
            colocated
        >
    > coloFaceFlux
    (
        const meshField<Type,colocated>&
    );

public:

    // Runtime type information
    TypeName("midPoint");

    // Construct from stream and mesh
    midPointFaceFluxScheme(const fvMesh&, Istream& = scheme::nullStream);

    // Face flux of a colocated scalar field
    virtual tmp<colocatedScalarFaceField> faceFlux
    (
        const colocatedVectorField&
    );

    // Face flux of a colocated tensor field
    virtual tmp<colocatedVectorFaceField> faceFlux
    (
        const colocatedTensorField&
    );

    // Face flux of a staggered scalar field
    virtual tmp<staggeredScalarFaceField> faceFlux
    (
        const staggeredScalarField&
    );
};

}

}

}

#ifdef NoRepository
    #include "midPointFaceFluxScheme.C"
#endif

#endif
