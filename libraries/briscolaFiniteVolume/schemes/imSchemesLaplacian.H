#ifndef imSchemesLaplacian_H
#define imSchemesLaplacian_H

#include "laplacianScheme.H"

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
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")",
            field.fvMsh()
        )->imLaplacian(lambdaPtr,field,factor);
    }

    // Coefficient as reference
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const meshField<lowerFaceScalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            "laplacian("+lambda.name()+","+field.name()+")",
            field.fvMsh()
        )->imLaplacian(lambda,field,factor);
    }

    // Coefficient as tmp
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const tmp<meshField<lowerFaceScalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            "laplacian("+tLambda->name()+","+field.name()+")",
            field.fvMsh()
        )->imLaplacian(tLambda,field,factor);
    }

    // Constant coefficient
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const scalar lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            "laplacian(const,"+field.name()+")",
            field.fvMsh()
        )->imLaplacian(lambda,field,factor);
    }

    // No coefficient
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            "laplacian("+field.name()+")",
            field.fvMsh()
        )->imLaplacian(field,factor);
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
    inline void laplacian
    (
        linearSystem<SType,Type,MeshType>& sys,
        const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")",
            field.fvMsh()
        )->template imLaplacian<OpType>(sys,lambdaPtr,field,factor);
    }

    // Coefficient as reference
    template
    <
        template<class> class OpType,
        class SType,
        class Type,
        class MeshType
    >
    inline void laplacian
    (
        linearSystem<SType,Type,MeshType>& sys,
        const meshField<lowerFaceScalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        laplacianScheme<SType,Type,MeshType>::New
        (
            "laplacian("+lambda.name()+","+field.name()+")",
            field.fvMsh()
        )->template imLaplacian<OpType>(sys,lambda,field,factor);
    }

    // Coefficient as tmp
    template
    <
        template<class> class OpType,
        class SType,
        class Type,
        class MeshType
    >
    inline void laplacian
    (
        linearSystem<SType,Type,MeshType>& sys,
        const tmp<meshField<lowerFaceScalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        laplacianScheme<SType,Type,MeshType>::New
        (
            "laplacian("+tLambda->name()+","+field.name()+")",
            field.fvMsh()
        )->template imLaplacian<OpType>(sys,tLambda,field,factor);
    }

    // Constant coefficient
    template
    <
        template<class> class OpType,
        class SType,
        class Type,
        class MeshType
    >
    inline void laplacian
    (
        linearSystem<SType,Type,MeshType>& sys,
        const scalar lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        laplacianScheme<SType,Type,MeshType>::New
        (
            "laplacian(const,"+field.name()+")",
            field.fvMsh()
        )->template imLaplacian<OpType>(sys,lambda,field,factor);
    }

    // No coefficient
    template
    <
        template<class> class OpType,
        class SType,
        class Type,
        class MeshType
    >
    inline void laplacian
    (
        linearSystem<SType,Type,MeshType>& sys,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            "laplacian("+field.name()+")",
            field.fvMsh()
        )->template imLaplacian<OpType>(sys,field,factor);
    }

// Default symmetric stencil type on colocated meshes

    // Return

        // Coefficient as pointer
        template<class Type>
        inline tmp<linearSystem<symmStencil,Type,colocated>> laplacian
        (
            const meshField<lowerFaceScalar,colocated>* lambdaPtr,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<symmStencil,Type,colocated>
            (
                lambdaPtr,
                field,
                factor
            );
        }

        // Coefficient as reference
        template<class Type>
        inline tmp<linearSystem<symmStencil,Type,colocated>> laplacian
        (
            const meshField<lowerFaceScalar,colocated>& lambda,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<symmStencil,Type,colocated>
            (
                lambda,
                field,
                factor
            );
        }

        // Coefficient as tmp
        template<class Type>
        inline tmp<linearSystem<symmStencil,Type,colocated>> laplacian
        (
            const tmp<meshField<lowerFaceScalar,colocated>>& tLambda,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<symmStencil,Type,colocated>
            (
                tLambda,
                field,
                factor
            );
        }

        // Constant coefficient
        template<class Type>
        inline tmp<linearSystem<symmStencil,Type,colocated>> laplacian
        (
            const scalar lambda,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<symmStencil,Type,colocated>
            (
                lambda,
                field,
                factor
            );
        }

        // No coefficient
        template<class Type>
        inline tmp<linearSystem<symmStencil,Type,colocated>> laplacian
        (
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<symmStencil,Type,colocated>(field,factor);
        }

// Default stencil type on staggered meshes

    // Return

        // Coefficient as pointer
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const meshField<lowerFaceScalar,staggered>* lambdaPtr,
            const meshField<Type,staggered>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil,Type,staggered>
            (
                lambdaPtr,
                field,
                factor
            );
        }

        // Coefficient as reference
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const meshField<lowerFaceScalar,staggered>& lambda,
            const meshField<Type,staggered>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil,Type,staggered>
            (
                lambda,
                field,
                factor
            );
        }

        // Coefficient as tmp
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const tmp<meshField<lowerFaceScalar,staggered>>& tLambda,
            const meshField<Type,staggered>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil,Type,staggered>
            (
                tLambda,
                field,
                factor
            );
        }

        // Constant coefficient
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const scalar lambda,
            const meshField<Type,staggered>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil,Type,staggered>
            (
                lambda,
                field,
                factor
            );
        }

        // No coefficient
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const meshField<Type,staggered>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil,Type,staggered>(field,factor);
        }
}

}

}

}

#endif
