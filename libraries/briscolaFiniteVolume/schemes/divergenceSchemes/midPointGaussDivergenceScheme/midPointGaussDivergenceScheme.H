#ifndef midPointGaussDivergenceScheme_H
#define midPointGaussDivergenceScheme_H

#include "divergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class midPointGaussDivergenceScheme
:
    public divergenceScheme<stencil,Type,MeshType>
{
protected:

    // Implicit divergence scheme
    template<template<class> class OpType>
    void imDiv
    (
        linearSystem<stencil,Type,MeshType>& sys,
        const meshField<lowerFaceScalar,MeshType>& phi,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    );

public:

    // Runtime type information
    TypeName("midPointGauss");

    // Construct from dictionary and mesh
    midPointGaussDivergenceScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    midPointGaussDivergenceScheme(const fvMesh&);

    #define DIVERGENCE(OPNAME)                                                 \
                                                                               \
    virtual void imDiv                                                         \
    (                                                                          \
        linearSystem<stencil,Type,MeshType>& sys,                              \
        const meshField<lowerFaceScalar,MeshType>& phi,                        \
        const meshField<Type,MeshType>& field,                                 \
        const OPNAME<Type>& op,                                                \
        const scalar factor = 1.0                                              \
    )                                                                          \
    {                                                                          \
        this->imDiv<OPNAME>(sys,phi,field,factor);                             \
    }

    DIVERGENCE(eqOp)
    DIVERGENCE(eqMinusOp)
    DIVERGENCE(plusEqOp)
    DIVERGENCE(minusEqOp)

    #undef DIVERGENCE

    virtual tmp<meshField<Type,MeshType>> exDiv
    (
        const meshField<lowerFaceScalar,MeshType>& phi,
        const meshField<Type,MeshType>& field
    );
};

}

}

}

#ifdef NoRepository
    #include "midPointGaussDivergenceScheme.C"
#endif

#endif
