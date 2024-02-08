#ifndef imSchemesDivergence_H
#define imSchemesDivergence_H

#include "divergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace im
{

// Return

    // Coefficient as reference
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> div
    (
        const meshField<lowerFaceScalar,MeshType>& phi,
        const meshField<Type,MeshType>& field
    )
    {
        return divergenceScheme<SType,Type,MeshType>::New
        (
            "div("+phi.name()+","+field.name()+")",
            field.fvMsh()
        )->imDiv(phi,field);
    }

    // Coefficient as tmp
    template<class SType, class Type, class MeshType>
    inline tmp<linearSystem<SType,Type,MeshType>> div
    (
        const tmp<meshField<lowerFaceScalar,MeshType>>& tPhi,
        const meshField<Type,MeshType>& field
    )
    {
        return divergenceScheme<SType,Type,MeshType>::New
        (
            "div("+tPhi->name()+","+field.name()+")",
            field.fvMsh()
        )->imDiv(tPhi,field);
    }

// Default stencil type on colocated and staggered meshes

    // Return

        // Coefficient as reference
        template<class Type, class MeshType>
        inline tmp<linearSystem<stencil,Type,MeshType>> div
        (
            const meshField<lowerFaceScalar,MeshType>& phi,
            const meshField<Type,MeshType>& field
        )
        {
            return div<stencil>(phi,field);
        }

        // Coefficient as tmp
        template<class Type, class MeshType>
        inline tmp<linearSystem<stencil,Type,MeshType>> div
        (
            const tmp<meshField<lowerFaceScalar,MeshType>>& tPhi,
            const meshField<Type,MeshType>& field
        )
        {
            return div<stencil>(tPhi,field);
        }
}

}

}

}

#endif
