#ifndef exSchemesDdt_H
#define exSchemesDdt_H

#include "ddtScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Return functions

    // Coefficient as pointer

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> ddt
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
        )->exDdt(lambdaPtr,field,factor);
    }

    // Coefficient as reference

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> ddt
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
        )->exDdt(lambda,field,factor);
    }

    // Coefficient as tmp

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> ddt
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
        )->exDdt(tLambda,field,factor);
    }

    // Constant coefficient

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> ddt
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
        )->exDdt(lambda,field,factor);
    }

    // No coefficient

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> ddt
    (
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            "ddt("+field.name()+")",
            field.fvMsh()
        )->exDdt(field,factor);
    }

}

}

}

}

#endif
