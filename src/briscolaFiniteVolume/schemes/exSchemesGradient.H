#ifndef exSchemesGradient_H
#define exSchemesGradient_H

#include "gradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Centered gradient

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const meshField<Type,MeshType>& field
)
{
    return gradientScheme<Type,MeshType>::New
    (
        field.fvMsh(),
        "grad("+field.name()+")"
    )->grad(field);
}

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const tmp<meshField<Type,MeshType>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    typedef typename outerProduct<vector,Type>::type GradType;

    tmp<meshField<GradType,MeshType>> tGrad = grad(tField());

    if (tField.isTmp())
        tField.clear();

    return tGrad;
}

// Centered gradient from a face field

using ::Foam::briscola::fv::grad;

// Staggered gradient from a colocated field

template<class Type>
inline tmp<meshField<Type,staggered>>
stagGrad
(
    const meshField<Type,colocated>& field
)
{
    return gradientScheme<Type,colocated>::New
    (
        field.fvMsh(),
        "stagGrad("+field.name()+")"
    )->stagGrad(field);
}

template<class Type>
inline tmp<meshField<Type,staggered>>
stagGrad
(
    const tmp<meshField<Type,colocated>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<Type,staggered>> tStagGrad = stagGrad(tField());

    if (tField.isTmp())
        tField.clear();

    return tStagGrad;
}

}

}

}

}

#endif
