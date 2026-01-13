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
            const faceField<scalar,MeshType>& lambda,
            const meshField<Type,MeshType>& field
        )
        {
            return this->imLaplacian(&lambda, field);
        }

        // Coefficient as tmp

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            const tmp<faceField<scalar,MeshType>>& tLambda,
            const meshField<Type,MeshType>& field
        )
        {
            tmp<linearSystem<SType,Type,MeshType>> tSys =
                this->imLaplacian(&tLambda(), field);

            if (tLambda.isTmp())
                tLambda.clear();

            return tSys;
        }

        // Constant coefficient

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            const scalar lambda,
            const meshField<Type,MeshType>& field
        )
        {
            return lambda*this->imLaplacian(nullptr, field);
        }

        // No coefficient

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::imLaplacian
        (
            const meshField<Type,MeshType>& field
        )
        {
            return this->imLaplacian(nullptr, field);
        }

// Explicit

    // Return

        // Coefficient as reference

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::exLaplacian
        (
            const faceField<scalar,MeshType>& lambda,
            const meshField<Type,MeshType>& field
        )
        {
            return this->exLaplacian(&lambda, field);
        }

        // Coefficient as tmp

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::exLaplacian
        (
            const tmp<faceField<scalar,MeshType>>& tLambda,
            const meshField<Type,MeshType>& field
        )
        {
            tmp<meshField<Type,MeshType>> tLap =
                this->exLaplacian(&tLambda(), field);

            if (tLambda.isTmp())
                tLambda.clear();

            return tLap;
        }

        // Constant coefficient

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::exLaplacian
        (
            const scalar lambda,
            const meshField<Type,MeshType>& field
        )
        {
            return lambda*this->exLaplacian(nullptr, field);
        }

        // No coefficient

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>>
        laplacianScheme<SType,Type,MeshType>::exLaplacian
        (
            const meshField<Type,MeshType>& field
        )
        {
            return this->exLaplacian(nullptr, field);
        }
}

}

}
