namespace Foam
{

namespace briscola
{

namespace fv
{

// Implicit

    // Return

        // Coefficient as pointer

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            tmp<linearSystem<SType,Type,MeshType>> tSys
            (
                new linearSystem<SType,Type,MeshType>
                (
                    const_cast<meshField<Type,MeshType>&>(field)
                )
            );

            linearSystem<SType,Type,MeshType>& sys = tSys.ref();

            this->imLaplacian
            (
                sys,
                lambdaPtr,
                field,
                eqOp<Type>(),
                factor
            );

            return tSys;
        }

        // Coefficient as reference

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            const meshField<lowerFaceScalar,MeshType>& lambda,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            return this->imLaplacian(&lambda, field, factor);
        }

        // Coefficient as tmp

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            const tmp<meshField<lowerFaceScalar,MeshType>>& tLambda,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            tmp<linearSystem<SType,Type,MeshType>> tSys
            (
                this->imLaplacian(&tLambda(), field, factor)
            );

            const_cast<tmp<meshField<lowerFaceScalar,MeshType>>&>(tLambda)
           .clear();

            return tSys;
        }

        // Constant coefficient

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            const scalar lambda,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            return this->imLaplacian(nullptr, field, lambda*factor);
        }

        // No coefficient

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            return this->imLaplacian(nullptr, field, factor);
        }

    // Manipulate linear system

        // Coefficient as pointer

        template<class SType, class Type, class MeshType>
        template<template<class> class OpType>
        inline void laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            linearSystem<SType,Type,MeshType>& sys,
            const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            this->imLaplacian(sys, lambdaPtr, field, OpType<Type>(), factor);
        }

        // Coefficient as reference

        template<class SType, class Type, class MeshType>
        template<template<class> class OpType>
        inline void laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            linearSystem<SType,Type,MeshType>& sys,
            const meshField<lowerFaceScalar,MeshType>& lambda,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            this->imLaplacian(sys, &lambda, field, OpType<Type>(), factor);
        }

        // Coefficient as tmp

        template<class SType, class Type, class MeshType>
        template<template<class> class OpType>
        inline void laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            linearSystem<SType,Type,MeshType>& sys,
            const tmp<meshField<lowerFaceScalar,MeshType>>& tLambda,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            this->imLaplacian(sys, &tLambda(), field, OpType<Type>(), factor);
            const_cast<tmp<meshField<lowerFaceScalar,MeshType>>&>(tLambda)
           .clear();
        }

        // Constant coefficient

        template<class SType, class Type, class MeshType>
        template<template<class> class OpType>
        inline void laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            linearSystem<SType,Type,MeshType>& sys,
            const scalar lambda,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            this->imLaplacian
            (
                sys,
                nullptr,
                field,
                OpType<Type>(),
                factor*lambda
            );
        }

        // No coefficient

        template<class SType, class Type, class MeshType>
        template<template<class> class OpType>
        inline void laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            linearSystem<SType,Type,MeshType>& sys,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            this->imLaplacian(sys, nullptr, field, OpType<Type>(), factor);
        }

// Explicit

    // Return

        // Coefficient as reference

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::exLaplacian
        (
            const meshField<lowerFaceScalar,MeshType>& lambda,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            return this->exLaplacian(&lambda, field, factor);
        }

        // Coefficient as tmp

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::exLaplacian
        (
            const tmp<meshField<lowerFaceScalar,MeshType>>& tLambda,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            tmp<meshField<Type,MeshType>> tLap
            (
                this->exLaplacian(&tLambda(), field, factor)
            );

            const_cast<tmp<meshField<lowerFaceScalar,MeshType>>&>(tLambda)
           .clear();

            return tLap;
        }

        // Constant coefficient

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::exLaplacian
        (
            const scalar lambda,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            return this->exLaplacian(nullptr, field, lambda*factor);
        }

        // No coefficient

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::exLaplacian
        (
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            return this->exLaplacian(nullptr, field, factor);
        }
}

}

}
