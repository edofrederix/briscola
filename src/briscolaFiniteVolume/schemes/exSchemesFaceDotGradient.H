#ifndef exSchemesFaceDotGradient_H
#define exSchemesFaceDotGradient_H

#include "faceDotGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Face-dot gradient

template<class Type, class MeshType>
inline tmp<faceField<Type,MeshType>> faceDotGrad
(
    const meshField<Type,MeshType>& field
)
{
    return faceDotGradientScheme<Type,MeshType>::New
    (
        field.fvMsh(),
        "faceDotGrad("+field.name()+")"
    )->faceDotGrad(field);
}

template<class Type, class MeshType>
inline tmp<faceField<Type,MeshType>> faceDotGrad
(
    const tmp<meshField<Type,MeshType>>& tField
)
{
    if (tField.isTmp())
        const_cast<tmp<meshField<Type,MeshType>>&>(tField)
            ->correctBoundaryConditions();

    tmp<faceField<Type,MeshType>> tFaceDotGrad =
        faceDotGrad(tField());

    if (tField.isTmp())
        tField.clear();

    return tFaceDotGrad;
}

}

}

}

}

#endif
