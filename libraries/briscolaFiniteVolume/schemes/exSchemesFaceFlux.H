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

}

}

}

}

#endif
