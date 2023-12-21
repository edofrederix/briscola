namespace Foam
{

namespace briscola
{

namespace fv
{

// Implicit

    // Return

        // Flux as reference

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        divergenceScheme<SType,Type,MeshType>::imDiv
        (
            const meshField<lowerFaceScalar,MeshType>& phi,
            const meshField<Type,MeshType>& field
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

            this->imDiv
            (
                sys,
                phi,
                field,
                eqOp<Type>()
            );

            return tSys;
        }

        // Flux as tmp

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        divergenceScheme<SType,Type,MeshType>::imDiv
        (
            const tmp<meshField<lowerFaceScalar,MeshType>>& tPhi,
            const meshField<Type,MeshType>& field
        )
        {
            tmp<linearSystem<SType,Type,MeshType>> tSys
            (
                this->imDiv(tPhi(), field)
            );

            const_cast<tmp<meshField<lowerFaceScalar,MeshType>>&>(tPhi)
           .clear();

            return tSys;
        }

    // Manipulate linear system

        // Flux as reference

        template<class SType, class Type, class MeshType>
        template<template<class> class OpType>
        inline void divergenceScheme<SType,Type,MeshType>::imDiv
        (
            linearSystem<SType,Type,MeshType>& sys,
            const meshField<lowerFaceScalar,MeshType>& phi,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            this->imDiv(sys, phi, field, OpType<Type>(), factor);
        }

        // Flux as tmp

        template<class SType, class Type, class MeshType>
        template<template<class> class OpType>
        inline void divergenceScheme<SType,Type,MeshType>::imDiv
        (
            linearSystem<SType,Type,MeshType>& sys,
            const tmp<meshField<lowerFaceScalar,MeshType>>& tPhi,
            const meshField<Type,MeshType>& field,
            const scalar factor
        )
        {
            this->imDiv(sys, tPhi(), field, OpType<Type>(), factor);
            const_cast<tmp<meshField<lowerFaceScalar,MeshType>>&>(tPhi)
           .clear();
        }

// Explicit

    // Return

        // Flux as tmp

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>>
        divergenceScheme<SType,Type,MeshType>::exDiv
        (
            const tmp<meshField<lowerFaceScalar,MeshType>>& tPhi,
            const meshField<Type,MeshType>& field
        )
        {
            tmp<meshField<Type,MeshType>> tLap
            (
                this->exDiv(tPhi(), field)
            );

            const_cast<tmp<meshField<lowerFaceScalar,MeshType>>&>(tPhi)
           .clear();

            return tLap;
        }
}

}

}
