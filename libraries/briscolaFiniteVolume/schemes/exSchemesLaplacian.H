#ifndef exSchemesLaplacian_H
#define exSchemesLaplacian_H

#include "laplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// General functions

    // Return

        // Coefficient as pointer

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>> laplacian
        (
            const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
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
            )->exLaplacian(lambdaPtr,field,factor);
        }

        // Coefficient as reference

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>> laplacian
        (
            const meshField<lowerFaceScalar,MeshType>& lambda,
            const meshField<Type,MeshType>& field,
            const scalar factor = 1.0
        )
        {
            return laplacianScheme<SType,Type,MeshType>::New
            (
                field.fvMsh(),
                "laplacian("+lambda.name()+","+field.name()+")"
            )->exLaplacian(lambda,field,factor);
        }

        // Coefficient as tmp

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>> laplacian
        (
            const tmp<meshField<lowerFaceScalar,MeshType>>& tLambda,
            const meshField<Type,MeshType>& field,
            const scalar factor = 1.0
        )
        {
            return laplacianScheme<SType,Type,MeshType>::New
            (
                field.fvMsh(),
                "laplacian("+tLambda->name()+","+field.name()+")"
            )->exLaplacian(tLambda,field,factor);
        }

        // Constant coefficient

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>> laplacian
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
            )->exLaplacian(lambda,field,factor);
        }

        // No coefficient

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>> laplacian
        (
            const meshField<Type,MeshType>& field,
            const scalar factor = 1.0
        )
        {
            return laplacianScheme<SType,Type,MeshType>::New
            (
                field.fvMsh(),
                "laplacian("+field.name()+")"
            )->exLaplacian(field,factor);
        }

// Default symmetric stencil type on colocated meshes

    // Return

        // Coefficient as pointer

        template<class Type>
        inline tmp<meshField<Type,colocated>> laplacian
        (
            const meshField<lowerFaceScalar,colocated>* lambdaPtr,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<symmStencil>(lambdaPtr,field,factor);
        }

        // Coefficient as reference

        template<class Type>
        inline tmp<meshField<Type,colocated>> laplacian
        (
            const meshField<lowerFaceScalar,colocated>& lambda,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<symmStencil>(lambda,field,factor);
        }

        // Coefficient as tmp

        template<class Type>
        inline tmp<meshField<Type,colocated>> laplacian
        (
            const tmp<meshField<lowerFaceScalar,colocated>>& tLambda,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<symmStencil>(tLambda,field,factor);
        }

        // Coefficient as reference

        template<class Type>
        inline tmp<meshField<Type,colocated>> laplacian
        (
            const scalar lambda,
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<symmStencil>(lambda,field,factor);
        }

        // No coefficient

        template<class Type>
        inline tmp<meshField<Type,colocated>> laplacian
        (
            const meshField<Type,colocated>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<symmStencil>(field,factor);
        }

// Default stencil type on staggered meshes

    // Return

        // Coefficient as pointer

        template<class Type>
        inline tmp<meshField<Type,staggered>> laplacian
        (
            const meshField<lowerFaceScalar,staggered>* lambdaPtr,
            const meshField<Type,staggered>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil>(lambdaPtr,field,factor);
        }

        // Coefficient as reference

        template<class Type>
        inline tmp<meshField<Type,staggered>> laplacian
        (
            const meshField<lowerFaceScalar,staggered>& lambda,
            const meshField<Type,staggered>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil>(lambda,field,factor);
        }

        // Coefficient as tmp

        template<class Type>
        inline tmp<meshField<Type,staggered>> laplacian
        (
            const tmp<meshField<lowerFaceScalar,staggered>>& tLambda,
            const meshField<Type,staggered>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil>(tLambda,field,factor);
        }

        // Constant coefficient

        template<class Type>
        inline tmp<meshField<Type,staggered>> laplacian
        (
            const scalar lambda,
            const meshField<Type,staggered>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil>(lambda,field,factor);
        }

        // No coefficient

        template<class Type>
        inline tmp<meshField<Type,staggered>> laplacian
        (
            const meshField<Type,staggered>& field,
            const scalar factor = 1.0
        )
        {
            return laplacian<stencil>(field,factor);
        }
}

}

}

}

#endif
