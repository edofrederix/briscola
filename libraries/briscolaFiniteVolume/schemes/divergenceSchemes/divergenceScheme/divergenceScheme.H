#ifndef divergenceScheme_H
#define divergenceScheme_H

#include "scheme.H"
#include "linearSystem.H"
#include "stencil.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class divergenceScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("divergenceScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        divergenceScheme,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    divergenceScheme(const dictionary&, const fvMesh&);

    // Copy constructor
    divergenceScheme(const divergenceScheme<SType,Type,MeshType>&);

    // Selector
    static autoPtr<divergenceScheme<SType,Type,MeshType>> New
    (
        const word,
        const fvMesh&
    );

    // Destructor
    virtual ~divergenceScheme();

    // Pure virtual functions, to be defined by derived schemes. Cannot use
    // templating here because they are virtual.

        #define DIVERGENCE(OPNAME)                                             \
                                                                               \
        virtual void imDiv                                                     \
        (                                                                      \
            linearSystem<SType,Type,MeshType>& sys,                            \
            const meshField<lowerFaceScalar,MeshType>& phi,                    \
            const meshField<Type,MeshType>& field,                             \
            const OPNAME<Type>& op,                                            \
            const scalar factor = 1.0                                          \
        ) = 0;

        DIVERGENCE(eqOp)
        DIVERGENCE(eqMinusOp)
        DIVERGENCE(plusEqOp)
        DIVERGENCE(minusEqOp)

        #undef DIVERGENCE

        virtual tmp<meshField<Type,MeshType>> exDiv
        (
            const meshField<lowerFaceScalar,MeshType>& phi,
            const meshField<Type,MeshType>& field
        ) = 0;

    // Support functions relying on the pure virtual functions above

        // Implicit

            // Return

                // Flux as reference
                inline tmp<linearSystem<SType,Type,MeshType>> imDiv
                (
                    const meshField<lowerFaceScalar,MeshType>& phi,
                    const meshField<Type,MeshType>& field
                );

                // Flux as tmp
                inline tmp<linearSystem<SType,Type,MeshType>> imDiv
                (
                    const tmp<meshField<lowerFaceScalar,MeshType>>& tPhi,
                    const meshField<Type,MeshType>& field
                );

            // Manipulate linear system

                // Flux as reference
                template<template<class> class OpType>
                inline void imDiv
                (
                    linearSystem<SType,Type,MeshType>& sys,
                    const meshField<lowerFaceScalar,MeshType>& phi,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // Flux as tmp
                template<template<class> class OpType>
                inline void imDiv
                (
                    linearSystem<SType,Type,MeshType>& sys,
                    const tmp<meshField<lowerFaceScalar,MeshType>>& tPhi,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

        // Explicit

            // Return

                // Flux as tmp
                inline tmp<meshField<Type,MeshType>> exDiv
                (
                    const tmp<meshField<lowerFaceScalar,MeshType>>& tPhi,
                    const meshField<Type,MeshType>& field
                );
};

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> explicitDiv
(
    const meshField<LowerFaceSpace<Type>,MeshType>&
);

template<class Type>
tmp<meshField<Type,colocated>> explicitColoDiv
(
    const meshField<Type,staggered>&
);

}

}

}

#include "divergenceSchemeI.H"

#ifdef NoRepository
    #include "divergenceScheme.C"
#endif

#endif
