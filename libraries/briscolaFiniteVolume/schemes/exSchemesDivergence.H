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

    using ::Foam::briscola::fv::div;

    // Colocated divergence from staggered field

    using ::Foam::briscola::fv::coloDiv;

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
                field.fvMsh(),
                "div("+phi.name()+","+field.name()+")"
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
                field.fvMsh(),
                "div("+tPhi->name()+","+field.name()+")"
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
