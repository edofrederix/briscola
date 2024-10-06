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
            field.fvMsh(),
            "div("+phi.name()+","+field.name()+")"
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
            field.fvMsh(),
            "div("+tPhi->name()+","+field.name()+")"
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
