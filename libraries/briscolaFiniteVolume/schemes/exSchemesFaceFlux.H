#ifndef exSchemesFaceFlux_H
#define exSchemesFaceFlux_H

#include "faceFluxScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Colocated face flux from staggered field

using ::Foam::briscola::fv::coloFaceFlux;

// Colocated face flux from vector

inline tmp<meshField<lowerFaceScalar,colocated>> faceFlux
(
    const colocatedVectorField& field
)
{
    return faceFluxScheme::New
    (
        "flux("+field.name()+")",
        field.fvMsh()
    )->faceFlux(field);
}

inline tmp<meshField<lowerFaceScalar,colocated>> faceFlux
(
    const tmp<colocatedVectorField>& tField
)
{
    tmp<meshField<lowerFaceScalar,colocated>> tFaceFlux
    (
        faceFlux(tField())
    );

    tField.clear();

    return tFaceFlux;
}

// Staggered face flux from scalar

inline tmp<meshField<lowerFaceScalar,staggered>> faceFlux
(
    const staggeredScalarField& field
)
{
    return faceFluxScheme::New
    (
        "flux("+field.name()+")",
        field.fvMsh()
    )->faceFlux(field);
}

inline tmp<meshField<lowerFaceScalar,staggered>> faceFlux
(
    const tmp<staggeredScalarField>& tField
)
{
    tmp<meshField<lowerFaceScalar,staggered>> tFaceFlux
    (
        faceFlux(tField())
    );

    tField.clear();

    return tFaceFlux;
}

}

}

}

}

#endif
