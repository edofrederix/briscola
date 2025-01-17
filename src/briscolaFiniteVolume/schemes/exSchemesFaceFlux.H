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

inline tmp<meshField<faceScalar,colocated>> faceFlux
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

inline tmp<meshField<faceScalar,colocated>> faceFlux
(
    const tmp<colocatedVectorField>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<faceScalar,colocated>> tFaceFlux
    (
        faceFlux(tField())
    );

    if (tField.isTmp())
        tField.clear();

    return tFaceFlux;
}

// Colocated face flux from tensor

inline tmp<meshField<faceVector,colocated>> faceFlux
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

inline tmp<meshField<faceVector,colocated>> faceFlux
(
    const tmp<colocatedTensorField>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<faceVector,colocated>> tFaceFlux
    (
        faceFlux(tField())
    );

    if (tField.isTmp())
        tField.clear();

    return tFaceFlux;
}

// Staggered face flux from scalar

inline tmp<meshField<faceScalar,staggered>> faceFlux
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

inline tmp<meshField<faceScalar,staggered>> faceFlux
(
    const tmp<staggeredScalarField>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<faceScalar,staggered>> tFaceFlux
    (
        faceFlux(tField())
    );

    if (tField.isTmp())
        tField.clear();

    return tFaceFlux;
}

}

}

}

}

#endif
