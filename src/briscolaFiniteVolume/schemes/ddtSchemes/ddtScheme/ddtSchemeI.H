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
        const meshField<Type,MeshType>& field
    )
    {
        return this->imDdt(&lambda, field);
    }

    // Coefficient as tmp

    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>>
    ddtScheme<Type,MeshType>::imDdt
    (
        const tmp<meshField<scalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field
    )
    {
        tmp<linearSystem<diagStencil,Type,MeshType>> tSys =
            this->imDdt(&tLambda(), field);

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
        const meshField<Type,MeshType>& field
    )
    {
        return lambda*this->imDdt(nullptr, field);
    }

    // No coefficient

    template<class Type, class MeshType>
    inline tmp<linearSystem<diagStencil,Type,MeshType>>
    ddtScheme<Type,MeshType>::imDdt
    (
        const meshField<Type,MeshType>& field
    )
    {
        return this->imDdt(nullptr, field);
    }

// Explicit

    // Coefficient as reference

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>>
    ddtScheme<Type,MeshType>::exDdt
    (
        const meshField<scalar,MeshType>& lambda,
        const meshField<Type,MeshType>& field
    )
    {
        return this->exDdt(&lambda, field);
    }

    // Coefficient as tmp

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>>
    ddtScheme<Type,MeshType>::exDdt
    (
        const tmp<meshField<scalar,MeshType>>& tLambda,
        const meshField<Type,MeshType>& field
    )
    {
        tmp<meshField<Type,MeshType>> tDdt =
            this->exDdt(&tLambda(), field);

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
        const meshField<Type,MeshType>& field
    )
    {
        return lambda*this->exDdt(nullptr, field);
    }

    // No coefficient

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>>
    ddtScheme<Type,MeshType>::exDdt
    (
        const meshField<Type,MeshType>& field
    )
    {
        return this->exDdt(nullptr, field);
    }
}

}

}
