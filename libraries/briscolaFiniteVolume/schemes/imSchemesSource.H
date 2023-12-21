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
                const_cast<meshField<Type,MeshType>&>(field)
            )
        );

        linearSystem<diagStencil,Type,MeshType>& sys = tSys.ref();

        const meshField<scalar,MeshType>& cv =
            field.fvMsh().template metrics<MeshType>().cellVolumes();

        if (lambdaPtr)
        {
            sys.A() = factor*cv*(*lambdaPtr);
        }
        else
        {
            sys.A() = factor*cv;
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

// Manipulate linear system

    // Coefficient as pointer

    template
    <
        template<class> class OpType,
        class SType,
        class Type,
        class MeshType
    >
    inline void source
    (
        linearSystem<SType,Type,MeshType>& sys,
        const meshField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        if (lambdaPtr)
            const_cast<meshField<scalar,MeshType>&>(*lambdaPtr).restrict();

        const meshField<scalar,MeshType>& cv =
            field.fvMsh().template metrics<MeshType>().cellVolumes();

        meshField<SType,MeshType>& A = sys.A();
        meshField<Type,MeshType>& b = sys.b();

        if (isEqOp<OpType>() || isEqMinusOp<OpType>())
        {
            b = Zero;
        }

        OpType<scalar> bop;

        if (lambdaPtr)
        {
            forAllCells(A, l, d, i, j, k)
                bop(A(l,d,i,j,k)[0], factor*cv(l,d,i,j,k)*(*lambdaPtr)(l,d,i,j,k));
        }
        else
        {
            forAllCells(A, l, d, i, j, k)
                bop(A(l,d,i,j,k)[0], factor*cv(l,d,i,j,k));
        }

        if (lambdaPtr)
            const_cast<meshField<scalar,MeshType>&>(*lambdaPtr).makeShallow();
    }

    // Coefficient as reference

    template
    <
        template<class> class OpType,
        class SType,
        class Type,
        class MeshType
    >
    inline void source
    (
        linearSystem<SType,Type,MeshType>& sys,
        const meshField<scalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        source<OpType>(sys, &lambda, field, factor);
    }

    // Constant coefficient

    template
    <
        template<class> class OpType,
        class SType,
        class Type,
        class MeshType
    >
    inline void source
    (
        linearSystem<SType,Type,MeshType>& sys,
        const scalar lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        source<OpType>(sys, nullptr, field, factor*lambda);
    }
}

}

}

}

#endif
