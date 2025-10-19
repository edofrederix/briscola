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

inline tmp<faceField<scalar,colocated>> faceFlux
(
    const colocatedVectorField& field
)
{
    return faceFluxScheme::New
    (
        field.fvMsh(),
        "flux("+field.name()+")"
    )->faceFlux(field);
}

inline tmp<faceField<scalar,colocated>> faceFlux
(
    const tmp<colocatedVectorField>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<faceField<scalar,colocated>> tFaceFlux = faceFlux(tField());

    if (tField.isTmp())
        tField.clear();

    return tFaceFlux;
}

// Colocated face flux from tensor

inline tmp<faceField<vector,colocated>> faceFlux
(
    const colocatedTensorField& field
)
{
    return faceFluxScheme::New
    (
        field.fvMsh(),
        "flux("+field.name()+")"
    )->faceFlux(field);
}

inline tmp<faceField<vector,colocated>> faceFlux
(
    const tmp<colocatedTensorField>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<faceField<vector,colocated>> tFaceFlux = faceFlux(tField());

    if (tField.isTmp())
        tField.clear();

    return tFaceFlux;
}

// Staggered face flux from scalar

inline tmp<faceField<scalar,staggered>> faceFlux
(
    const staggeredScalarField& field
)
{
    return faceFluxScheme::New
    (
        field.fvMsh(),
        "flux("+field.name()+")"
    )->faceFlux(field);
}

inline tmp<faceField<scalar,staggered>> faceFlux
(
    const tmp<staggeredScalarField>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<faceField<scalar,staggered>> tFaceFlux = faceFlux(tField());

    if (tField.isTmp())
        tField.clear();

    return tFaceFlux;
}

}

}

}

}

#endif
