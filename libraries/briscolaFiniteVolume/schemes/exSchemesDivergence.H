#ifndef exSchemesDivergence_H
#define exSchemesDivergence_H

#include "divergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

// Functions not relying on scheme selection

    // Divergence term from flux

    template<class Type, class MeshType>
    inline tmp<meshField<Type,MeshType>> div
    (
        const meshField<LowerFaceSpace<Type>,MeshType>& phi
    )
    {
        return explicitDiv(phi);
    }

    // Colocated divergence from staggered field

    template<class Type>
    inline tmp<meshField<Type,colocated>> coloDiv
    (
        const meshField<Type,staggered>& field
    )
    {
        return explicitColoDiv(field);
    }

// General functions

    // Return

        // Coefficient as reference

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>> div
        (
            const meshField<lowerFaceScalar,MeshType>& phi,
            const meshField<Type,MeshType>& field
        )
        {
            return divergenceScheme<SType,Type,MeshType>::New
            (
                "div("+phi.name()+","+field.name()+")",
                field.fvMsh()
            )->exDiv(phi,field);
        }

        // Coefficient as tmp

        template<class SType, class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>> div
        (
            const tmp<meshField<lowerFaceScalar,MeshType>>& tPhi,
            const meshField<Type,MeshType>& field
        )
        {
            return divergenceScheme<SType,Type,MeshType>::New
            (
                "div("+tPhi->name()+","+field.name()+")",
                field.fvMsh()
            )->exDiv(tPhi,field);
        }

// Default stencil type on colocated meshes

    // Return

        // Coefficient as reference

        template<class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>> div
        (
            const meshField<lowerFaceScalar,MeshType>& phi,
            const meshField<Type,MeshType>& field
        )
        {
            return div<stencil>(phi,field);
        }

        // Coefficient as tmp

        template<class Type, class MeshType>
        inline tmp<meshField<Type,MeshType>> div
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
