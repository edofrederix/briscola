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
        const meshField<Type,MeshType>& field
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            field.fvMsh(),
            lambdaPtr
          ? "ddt("+lambdaPtr->name()+","+field.name()+")"
          : "ddt("+field.name()+")"
        )->exDdt(lambdaPtr,field);
    }

    // Coefficient as reference

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> ddt
    (
        const meshField<scalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            field.fvMsh(),
            "ddt("+lambda.name()+","+field.name()+")"
        )->exDdt(lambda,field);
    }

    // Coefficient as tmp

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> ddt
    (
        const tmp<meshField<scalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            field.fvMsh(),
            "ddt("+tLambda->name()+","+field.name()+")"
        )->exDdt(tLambda,field);
    }

    // Constant coefficient

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> ddt
    (
        const scalar lambda,
        const meshField<Type,MeshType>& field
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            field.fvMsh(),
            "ddt(const,"+field.name()+")"
        )->exDdt(lambda,field);
    }

    // No coefficient

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> ddt
    (
        const meshField<Type,MeshType>& field
    )
    {
        return ddtScheme<Type,MeshType>::New
        (
            field.fvMsh(),
            "ddt("+field.name()+")"
        )->exDdt(field);
    }

}

}

}

}

#endif
