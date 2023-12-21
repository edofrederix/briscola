#ifndef laplacianScheme_H
#define laplacianScheme_H

#include "scheme.H"
#include "linearSystem.H"
#include "stencil.H"
#include "ops.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class laplacianScheme
:
    public scheme
{
public:

    // Runtime type information
    TypeName("laplacianScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        laplacianScheme,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    laplacianScheme(const dictionary&, const fvMesh&);

    // Copy constructor
    laplacianScheme(const laplacianScheme<SType,Type,MeshType>&);

    // Selector
    static autoPtr<laplacianScheme<SType,Type,MeshType>> New
    (
        const word,
        const fvMesh&
    );

    // Destructor
    virtual ~laplacianScheme();

    // Pure virtual functions, to be defined by derived schemes. Cannot use
    // templating here because they are virtual.

        #define LAPLACIAN(OPNAME)                                              \
                                                                               \
        virtual void imLaplacian                                               \
        (                                                                      \
            linearSystem<SType,Type,MeshType>& sys,                            \
            const meshField<lowerFaceScalar,MeshType>* lambdaPtr,              \
            const meshField<Type,MeshType>& field,                             \
            const OPNAME<Type>& op,                                            \
            const scalar factor = 1.0                                          \
        ) = 0;

        LAPLACIAN(eqOp)
        LAPLACIAN(eqMinusOp)
        LAPLACIAN(plusEqOp)
        LAPLACIAN(minusEqOp)

        #undef LAPLACIAN

        virtual tmp<meshField<Type,MeshType>> exLaplacian
        (
            const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
            const meshField<Type,MeshType>& field,
            const scalar factor = 1.0
        ) = 0;

    // Support functions relying on the pure virtual functions above

        // Implicit

            // Return

                // Coefficient as pointer
                inline tmp<linearSystem<SType,Type,MeshType>> imLaplacian
                (
                    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // Coefficient as reference
                inline tmp<linearSystem<SType,Type,MeshType>> imLaplacian
                (
                    const meshField<lowerFaceScalar,MeshType>& lambda,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // Coefficient as tmp
                inline tmp<linearSystem<SType,Type,MeshType>> imLaplacian
                (
                    const tmp<meshField<lowerFaceScalar,MeshType>>& tLambda,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // Constant coefficient
                inline tmp<linearSystem<SType,Type,MeshType>> imLaplacian
                (
                    const scalar lambda,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // No coefficient
                inline tmp<linearSystem<SType,Type,MeshType>> imLaplacian
                (
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

            // Manipulate linear system

                // Coefficient as pointer
                template<template<class> class OpType>
                inline void imLaplacian
                (
                    linearSystem<SType,Type,MeshType>& sys,
                    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // Coefficient as reference
                template<template<class> class OpType>
                inline void imLaplacian
                (
                    linearSystem<SType,Type,MeshType>& sys,
                    const meshField<lowerFaceScalar,MeshType>& lambda,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // Coefficient as tmp
                template<template<class> class OpType>
                inline void imLaplacian
                (
                    linearSystem<SType,Type,MeshType>& sys,
                    const tmp<meshField<lowerFaceScalar,MeshType>>& tLambda,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // Constant coefficient
                template<template<class> class OpType>
                inline void imLaplacian
                (
                    linearSystem<SType,Type,MeshType>& sys,
                    const scalar lambda,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // No coefficient
                template<template<class> class OpType>
                inline void imLaplacian
                (
                    linearSystem<SType,Type,MeshType>& sys,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

        // Explicit

            // Return

                // Coefficient as reference
                inline tmp<meshField<Type,MeshType>> exLaplacian
                (
                    const meshField<lowerFaceScalar,MeshType>& lambda,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // Coefficient as tmp
                inline tmp<meshField<Type,MeshType>> exLaplacian
                (
                    const tmp<meshField<lowerFaceScalar,MeshType>>& tLambda,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // Constant coefficient
                inline tmp<meshField<Type,MeshType>> exLaplacian
                (
                    const scalar lambda,
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );

                // No coefficient
                inline tmp<meshField<Type,MeshType>> exLaplacian
                (
                    const meshField<Type,MeshType>& field,
                    const scalar factor = 1.0
                );
};

}

}

}

#include "laplacianSchemeI.H"

#ifdef NoRepository
    #include "laplacianScheme.C"
#endif

#endif
