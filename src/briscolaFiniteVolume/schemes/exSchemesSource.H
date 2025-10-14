#ifndef exSchemesSource_H
#define exSchemesSource_H

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Return

    // Coefficient as pointer

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> source
    (
        const meshField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field
    )
    {
        tmp<meshField<Type,MeshType>> tS =
            meshField<Type,MeshType>::New
            (
                lambdaPtr
              ? "source("+lambdaPtr->name()+","+field.name()+")"
              : "source("+field.name()+")",
                field.fvMsh()
            );

        meshField<Type,MeshType>& S = tS.ref();

        if (lambdaPtr)
        {
            S = (*lambdaPtr)*field;
        }
        else
        {
            S = field;
        }

        return tS;
    }

    // Coefficient as reference

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> source
    (
        const meshField<scalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field
    )
    {
        return source(&lambda, field);
    }

    // Coefficient as tmp

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> source
    (
        const tmp<meshField<scalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field
    )
    {
        tmp<meshField<Type,MeshType>> tS = source(tLambda(), field);

        if (tLambda.isTmp())
            tLambda.clear();

        return tS;
    }

    // Constant coefficient

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> source
    (
        const scalar lambda,
        const meshField<Type,MeshType>& field
    )
    {
        return lambda*source(nullptr, field);
    }
}

}

}

}

#endif
