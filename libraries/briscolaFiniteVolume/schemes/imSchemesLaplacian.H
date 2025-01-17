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
        const meshField<faceScalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            field.fvMsh(),
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")"
        )->imLaplacian(lambdaPtr,field,factor);
    }

    // Coefficient as reference
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const meshField<faceScalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            field.fvMsh(),
            "laplacian("+lambda.name()+","+field.name()+")"
        )->imLaplacian(lambda,field,factor);
    }

    // Coefficient as tmp
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const tmp<meshField<faceScalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            field.fvMsh(),
            "laplacian("+tLambda->name()+","+field.name()+")"
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
            field.fvMsh(),
            "laplacian(const,"+field.name()+")"
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
            field.fvMsh(),
            "laplacian("+field.name()+")"
        )->imLaplacian(field,factor);
    }

// Default stencil type on colocated meshes

    // Return

        // Coefficient as pointer
        template<class Type>
        inline tmp<linearSystem<stencil,Type,colocated>> laplacian
        (
            const meshField<faceScalar,colocated>* lambdaPtr,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil,Type,colocated>
            (
                lambdaPtr,
                field,
                factor
            );
        }

        // Coefficient as reference
        template<class Type>
        inline tmp<linearSystem<stencil,Type,colocated>> laplacian
        (
            const meshField<faceScalar,colocated>& lambda,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil,Type,colocated>
            (
                lambda,
                field,
                factor
            );
        }

        // Coefficient as tmp
        template<class Type>
        inline tmp<linearSystem<stencil,Type,colocated>> laplacian
        (
            const tmp<meshField<faceScalar,colocated>>& tLambda,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil,Type,colocated>
            (
                tLambda,
                field,
                factor
            );
        }

        // Constant coefficient
        template<class Type>
        inline tmp<linearSystem<stencil,Type,colocated>> laplacian
        (
            const scalar lambda,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil,Type,colocated>
            (
                lambda,
                field,
                factor
            );
        }

        // No coefficient
        template<class Type>
        inline tmp<linearSystem<stencil,Type,colocated>> laplacian
        (
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil,Type,colocated>(field,factor);
        }

// Default stencil type on staggered meshes

    // Return

        // Coefficient as pointer
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const meshField<faceScalar,staggered>* lambdaPtr,
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
            const meshField<faceScalar,staggered>& lambda,
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
            const tmp<meshField<faceScalar,staggered>>& tLambda,
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
