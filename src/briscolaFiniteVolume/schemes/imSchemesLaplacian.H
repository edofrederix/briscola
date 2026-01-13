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
        const faceField<scalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            field.fvMsh(),
            lambdaPtr
          ? "laplacian("+lambdaPtr->name()+","+field.name()+")"
          : "laplacian("+field.name()+")"
        )->imLaplacian(lambdaPtr,field);
    }

    // Coefficient as reference
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const faceField<scalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            field.fvMsh(),
            "laplacian("+lambda.name()+","+field.name()+")"
        )->imLaplacian(lambda,field);
    }

    // Coefficient as tmp
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const tmp<faceField<scalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            field.fvMsh(),
            "laplacian("+tLambda->name()+","+field.name()+")"
        )->imLaplacian(tLambda,field);
    }

    // Constant coefficient
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const scalar lambda,
        const meshField<Type,MeshType>& field
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            field.fvMsh(),
            "laplacian(const,"+field.name()+")"
        )->imLaplacian(lambda,field);
    }

    // No coefficient
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> laplacian
    (
        const meshField<Type,MeshType>& field
    )
    {
        return laplacianScheme<SType,Type,MeshType>::New
        (
            field.fvMsh(),
            "laplacian("+field.name()+")"
        )->imLaplacian(field);
    }

// Default stencil type on colocated meshes

    // Return

        // Coefficient as pointer
        template<class Type>
        inline tmp<linearSystem<stencil,Type,colocated>> laplacian
        (
            const faceField<scalar,colocated>* lambdaPtr,
            const meshField<Type,colocated>& field
        )
        {
            return laplacian<stencil,Type,colocated>
            (
                lambdaPtr,
                field
            );
        }

        // Coefficient as reference
        template<class Type>
        inline tmp<linearSystem<stencil,Type,colocated>> laplacian
        (
            const faceField<scalar,colocated>& lambda,
            const meshField<Type,colocated>& field
        )
        {
            return laplacian<stencil,Type,colocated>
            (
                lambda,
                field
            );
        }

        // Coefficient as tmp
        template<class Type>
        inline tmp<linearSystem<stencil,Type,colocated>> laplacian
        (
            const tmp<faceField<scalar,colocated>>& tLambda,
            const meshField<Type,colocated>& field
        )
        {
            return laplacian<stencil,Type,colocated>
            (
                tLambda,
                field
            );
        }

        // Constant coefficient
        template<class Type>
        inline tmp<linearSystem<stencil,Type,colocated>> laplacian
        (
            const scalar lambda,
            const meshField<Type,colocated>& field
        )
        {
            return lambda*laplacian<stencil,Type,colocated>(field);
        }

        // No coefficient
        template<class Type>
        inline tmp<linearSystem<stencil,Type,colocated>> laplacian
        (
            const meshField<Type,colocated>& field
        )
        {
            return laplacian<stencil,Type,colocated>(field);
        }

// Default stencil type on staggered meshes

    // Return

        // Coefficient as pointer
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const faceField<scalar,staggered>* lambdaPtr,
            const meshField<Type,staggered>& field
        )
        {
            return laplacian<stencil,Type,staggered>
            (
                lambdaPtr,
                field
            );
        }

        // Coefficient as reference
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const faceField<scalar,staggered>& lambda,
            const meshField<Type,staggered>& field
        )
        {
            return laplacian<stencil,Type,staggered>
            (
                lambda,
                field
            );
        }

        // Coefficient as tmp
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const tmp<faceField<scalar,staggered>>& tLambda,
            const meshField<Type,staggered>& field
        )
        {
            return laplacian<stencil,Type,staggered>
            (
                tLambda,
                field
            );
        }

        // Constant coefficient
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const scalar lambda,
            const meshField<Type,staggered>& field
        )
        {
            return lambda*laplacian<stencil,Type,staggered>(field);
        }

        // No coefficient
        template<class Type>
        inline tmp<linearSystem<stencil,Type,staggered>> laplacian
        (
            const meshField<Type,staggered>& field
        )
        {
            return laplacian<stencil,Type,staggered>(field);
        }
}

}

}

}

#endif
