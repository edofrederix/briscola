namespace Foam
{

namespace briscola
{

namespace fv
{

// Implicit

    // Coefficient as reference

    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>>
    ddtScheme<Type,MeshType>::imDdt
    (
        const meshField<scalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor
    )
    {
        return this->imDdt(&lambda, field, factor);
    }

    // Coefficient as tmp

    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>>
    ddtScheme<Type,MeshType>::imDdt
    (
        const tmp<meshField<scalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field,
        const scalar factor
    )
    {
        tmp<linearSystem<diagStencil,Type,MeshType>> tSys
        (
            this->imDdt(&tLambda(), field, factor)
        );

        if (tLambda.isTmp())
            tLambda.clear();

        return tSys;
    }

    // Constant coefficient

    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>>
    ddtScheme<Type,MeshType>::imDdt
    (
        const scalar lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor
    )
    {
        return this->imDdt(nullptr, field, lambda*factor);
    }

    // No coefficient

    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>>
    ddtScheme<Type,MeshType>::imDdt
    (
        const meshField<Type,MeshType>& field,
        const scalar factor
    )
    {
        return this->imDdt(nullptr, field, factor);
    }

// Explicit

    // Coefficient as reference

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>>
    ddtScheme<Type,MeshType>::exDdt
    (
        const meshField<scalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor
    )
    {
        return this->exDdt(&lambda, field, factor);
    }

    // Coefficient as tmp

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>>
    ddtScheme<Type,MeshType>::exDdt
    (
        const tmp<meshField<scalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field,
        const scalar factor
    )
    {
        tmp<meshField<Type,MeshType>> tDdt
        (
            this->exDdt(&tLambda(), field, factor)
        );

        if (tLambda.isTmp())
            tLambda.clear();

        return tDdt;
    }

    // Constant coefficient

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>>
    ddtScheme<Type,MeshType>::exDdt
    (
        const scalar lambda,
        const meshField<Type,MeshType>& field,
        const scalar factor
    )
    {
        return this->exDdt(nullptr, field, lambda*factor);
    }

    // No coefficient

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>>
    ddtScheme<Type,MeshType>::exDdt
    (
        const meshField<Type,MeshType>& field,
        const scalar factor
    )
    {
        return this->exDdt(nullptr, field, factor);
    }
}

}

}
