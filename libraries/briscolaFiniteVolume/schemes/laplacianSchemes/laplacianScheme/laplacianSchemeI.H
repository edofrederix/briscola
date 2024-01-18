namespace Foam
{

namespace briscola
{

namespace fv
{

// Implicit

    // Return

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
