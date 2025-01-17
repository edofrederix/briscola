namespace Foam
{

namespace briscola
{

namespace fv
{

// Implicit

    // Return

        // Flux as tmp

        template<class SType, class Type, class MeshType>
        inline tmp<linearSystem<SType,Type,MeshType>>
        divergenceScheme<SType,Type,MeshType>::imDiv
        (
            const tmp<meshField<faceScalar,MeshType>>& tPhi,
            const meshField<Type,MeshType>& field
        )
        {
            tmp<linearSystem<SType,Type,MeshType>> tSys
            (
                this->imDiv(tPhi(), field)
            );

            if (tPhi.isTmp())
                tPhi.clear();

            return tSys;
        }

// Explicit

    // Return

        // Flux as tmp

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>>
        divergenceScheme<SType,Type,MeshType>::exDiv
        (
            const tmp<meshField<faceScalar,MeshType>>& tPhi,
            const meshField<Type,MeshType>& field
        )
        {
            tmp<meshField<Type,MeshType>> tLap
            (
                this->exDiv(tPhi(), field)
            );

            if (tPhi.isTmp())
                tPhi.clear();

            return tLap;
        }
}

}

}
