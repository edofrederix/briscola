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
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
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
            S = factor*(*lambdaPtr)*field;
        }
        else
        {
            S = factor*field;
        }

        return tS;
    }

    // Coefficient as reference

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> source
    (
        const meshField<scalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return source(&lambda, field, factor);
    }

    // Coefficient as tmp

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> source
    (
        const tmp<meshField<scalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        tmp<meshField<Type,MeshType>> tS = source(tLambda(), field, factor);

        if (tLambda.isTmp())
            tLambda.clear();

        return tS;
    }

    // Constant coefficient

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> source
    (
        const scalar lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return source(nullptr, field, factor*lambda);
    }
}

}

}

}

#endif
