#ifndef symmStencilLinearGaussLaplacianScheme_H
#define symmStencilLinearGaussLaplacianScheme_H

#include "linearGaussLaplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
class symmStencilLinearGaussLaplacianScheme
:
    public linearGaussLaplacianScheme<symmStencil,Type,MeshType>
{
protected:

    // Implicit laplacian scheme
    template<template<class> class OpType>
    void imLaplacian
    (
        linearSystem<symmStencil,Type,MeshType>& sys,
        const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
        const meshField<Type,MeshType>& field,
        const scalar factor = 1.0
    );

public:

    // Construct from dictionary and mesh
    symmStencilLinearGaussLaplacianScheme(const dictionary&, const fvMesh&);

    // Construct from mesh only, with empty dictionary
    symmStencilLinearGaussLaplacianScheme(const fvMesh&);

    #define LAPLACIAN(OPNAME)                                                  \
                                                                               \
    virtual void imLaplacian                                                   \
    (                                                                          \
        linearSystem<symmStencil,Type,MeshType>& sys,                          \
        const meshField<lowerFaceScalar,MeshType>* lambdaPtr,                  \
        const meshField<Type,MeshType>& field,                                 \
        const OPNAME<Type>& op,                                                \
        const scalar factor = 1.0                                              \
    )                                                                          \
    {                                                                          \
        this->imLaplacian<OPNAME>(sys,lambdaPtr,field,factor);                 \
    }

    LAPLACIAN(eqOp)
    LAPLACIAN(eqMinusOp)
    LAPLACIAN(plusEqOp)
    LAPLACIAN(minusEqOp)

    #undef LAPLACIAN
};

}

}

}

#ifdef NoRepository
    #include "symmStencilLinearGaussLaplacianScheme.C"
#endif

#endif
