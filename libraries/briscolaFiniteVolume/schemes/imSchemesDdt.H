#ifndef imSchemesDdt_H
#define imSchemesDdt_H

#include "ddtScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace im
{

// Return functions

    // Coefficient as pointer
    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>> ddt
    (
        const meshField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            lambdaPtr
          ? "ddt("+lambdaPtr->name()+","+field.name()+")"
          : "ddt("+field.name()+")",
            field.fvMsh()
        )->imDdt(lambdaPtr,field,factor);
    }

    // Coefficient as reference
    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>> ddt
    (
        const meshField<scalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            "ddt("+lambda.name()+","+field.name()+")",
            field.fvMsh()
        )->imDdt(lambda,field,factor);
    }

    // Coefficient as tmp
    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>> ddt
    (
        const tmp<meshField<scalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            "ddt("+tLambda->name()+","+field.name()+")",
            field.fvMsh()
        )->imDdt(tLambda,field,factor);
    }

    // Constant coefficient
    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>> ddt
    (
        const scalar lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            "ddt(const,"+field.name()+")",
            field.fvMsh()
        )->imDdt(lambda,field,factor);
    }

    // No coefficient
    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>> ddt
    (
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            "ddt("+field.name()+")",
            field.fvMsh()
        )->imDdt(field,factor);
    }
}

}

}

}

#endif
