#ifndef exSchemes_H
#define exSchemes_H

#include "ddtScheme.H"
#include "divergenceScheme.H"
#include "faceFluxScheme.H"
#include "faceGradientScheme.H"
#include "gradientScheme.H"
#include "interpolationScheme.H"
#include "laplacianScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace ex
{

template<class Type, class MeshType>
inline tmp<meshField<Type,MeshType>> div
(
    const meshField<LowerFaceSpace<Type>,MeshType>& phi
)
{
    return explicitDiv(phi);
}

template<class Type, class MeshType>
inline tmp<meshField<Type,MeshType>> div
(
    const meshField<lowerFaceScalar,MeshType>& phi,
    meshField<Type,MeshType>& field
)
{
    return divergenceScheme<Type,MeshType>::New
    (
        "div("+phi.name()+","+field.name()+")",
        field.fvMsh()
    )->exDiv(phi,field);
}

template<class Type, class MeshType>
inline tmp<meshField<Type,MeshType>> laplacian
(
    const meshField<lowerFaceScalar,MeshType>& lambda,
    meshField<Type,MeshType>& field
)
{
    return laplacianScheme<Type,MeshType>::New
    (
        "laplacian("+lambda.name()+","+field.name()+")",
        field.fvMsh()
    )->exLaplacian(lambda,field);
}

template<class Type, class MeshType>
inline tmp<meshField<Type,MeshType>> laplacian
(
    const meshField<lowerFaceScalar,MeshType>* lambdaPtr,
    meshField<Type,MeshType>& field
)
{
    return laplacianScheme<Type,MeshType>::New
    (
        lambdaPtr
      ? "laplacian("+lambdaPtr->name()+","+ field.name()+")"
      : "laplacian("+field.name()+")",
        field.fvMsh()
    )->exLaplacian(lambdaPtr,field);
}

template<class Type, class MeshType>
inline tmp<meshField<Type,MeshType>> laplacian
(
    meshField<Type,MeshType>& field
)
{
    return laplacianScheme<Type,MeshType>::New
    (
        "laplacian("+field.name()+")",
        field.fvMsh()
    )->exLaplacian(field);
}

template<class Type>
inline tmp<meshField<Type,colocated>> coloDiv
(
    const meshField<Type,staggered>& field
)
{
    return explicitColoDiv(field);
}

inline tmp<meshField<lowerFaceScalar,colocated>> faceFlux
(
    const colocatedVectorField& field
)
{
    return faceFluxScheme::New
    (
        "flux("+field.name()+")",
        field.fvMsh()
    )->faceFlux(field);
}

inline tmp<meshField<lowerFaceScalar,staggered>> faceFlux
(
    const staggeredScalarField& field
)
{
    return faceFluxScheme::New
    (
        "flux("+field.name()+")",
        field.fvMsh()
    )->faceFlux(field);
}

template<class Type, class MeshType>
inline tmp<meshField<LowerFaceSpace<Type>,MeshType>> faceGrad
(
    const meshField<Type,MeshType>& field
)
{
    return faceGradientScheme<Type,MeshType>::New
    (
        "faceGrad("+field.name()+")",
        field.fvMsh()
    )->faceGrad(field);
}

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const meshField<Type,MeshType>& field
)
{
    return gradientScheme<Type,MeshType>::New
    (
        "grad("+field.name()+")",
        field.fvMsh()
    )->grad(field);
}

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const meshField<FaceSpace<Type>,MeshType>& field
)
{
    return centerGrad(field);
}

template<class Type>
inline tmp<meshField<Type,staggered>>
stagGrad
(
    const meshField<Type,colocated>& field
)
{
    return gradientScheme<Type,colocated>::New
    (
        "stagGrad("+field.name()+")",
        field.fvMsh()
    )->stagGrad(field);
}

template<class Type, class MeshType>
inline tmp<meshField<LowerFaceSpace<Type>,MeshType>> interp
(
    const meshField<Type,MeshType>& field
)
{
    return interpolationScheme<Type,MeshType>::New
    (
        "interp("+field.name()+")",
        field.fvMsh()
    )->interp(field);
}

template<class Type>
inline tmp<meshField<typename outerProduct<vector,Type>::type,colocated>>
reconstruct
(
    const meshField<Type,staggered>& field
)
{
    typedef typename outerProduct<vector,Type>::type ReconType;

    tmp<meshField<ReconType,colocated>> tRecon
    (
        new meshField<ReconType,colocated>
        (
            "reconstruct("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<ReconType,colocated>& Recon = tRecon.ref();

    Recon = Zero;

    const meshField<faceVector,colocated>& fn =
        field.fvMsh().template metrics<colocated>().faceNormals();

    forAllCells(Recon, i, j, k)
        Recon(i,j,k) =
            0.5
          * (
              - field(0,i,  j,  k  ) * fn(i,j,k).left()
              + field(0,i+1,j,  k  ) * fn(i,j,k).right()
              - field(1,i,  j,  k  ) * fn(i,j,k).bottom()
              + field(1,i,  j+1,k  ) * fn(i,j,k).top()
              - field(2,i,  j,  k  ) * fn(i,j,k).aft()
              + field(2,i,  j,  k+1) * fn(i,j,k).fore()
            );

    return tRecon;
}

}

}

}

}

#endif
