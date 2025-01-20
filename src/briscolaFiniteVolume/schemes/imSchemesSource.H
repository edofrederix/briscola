#ifndef imSchemesSource_H
#define imSchemesSource_H

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace im
{

// Return

    // Coefficient as pointer

    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>> source
    (
        const meshField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        if (lambdaPtr)
            const_cast<meshField<scalar,MeshType>&>(*lambdaPtr).restrict();

        tmp<linearSystem<diagStencil,Type,MeshType>> tSys
        (
            new linearSystem<diagStencil,Type,MeshType>
            (
                lambdaPtr
              ? word("source("+lambdaPtr->name()+","+field.name()+")")
              : word("source("+field.name()+")"),
                const_cast<meshField<Type,MeshType>&>(field)
            )
        );

        linearSystem<diagStencil,Type,MeshType>& sys = tSys.ref();

        const meshField<scalar,MeshType>& cv =
            field.fvMsh().template metrics<MeshType>().cellVolumes();

        sys.A() = factor*cv;

        if (lambdaPtr)
        {
            // The restricted lambda field has undefined ghosts. Manually
            // implement operation to avoid segmentation errors.

            forAllCells(sys.A(), l, d, i, j, k)
                sys.A()[l][d](i,j,k) *= (*lambdaPtr)[l][d](i,j,k);
        }

        sys.b() = Zero;

        if (lambdaPtr)
            const_cast<meshField<scalar,MeshType>&>(*lambdaPtr).makeShallow();

        return tSys;
    }

    // Coefficient as reference

    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>> source
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
    inline tmp<linearSystem<diagStencil,Type,MeshType>> source
    (
        const tmp<meshField<scalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        tmp<linearSystem<diagStencil,Type,MeshType>> tSource
        (
            source(tLambda(), field, factor)
        );

        if (tLambda.isTmp())
            tLambda.clear();

        return tSource;
    }

    // Constant coefficient

    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>> source
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
