#include "EulerDdtScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
EulerDdtScheme<Type,MeshType>::EulerDdtScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    ddtScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
EulerDdtScheme<Type,MeshType>::EulerDdtScheme(const fvMesh& fvMsh)
:
    ddtScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<linearSystem<diagStencil,Type,MeshType>>
EulerDdtScheme<Type,MeshType>::imDdt
(
    const meshField<scalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    if (lambdaPtr)
        const_cast<meshField<scalar,MeshType>&>(*lambdaPtr)
       .restrict();

    tmp<linearSystem<diagStencil,Type,MeshType>> tSys
    (
        new linearSystem<diagStencil,Type,MeshType>
        (
            const_cast<meshField<Type,MeshType>&>(field)
        )
    );

    linearSystem<diagStencil,Type,MeshType>& Sys = tSys.ref();

    const meshField<scalar,MeshType>& cv =
        this->fvMsh().template metrics<MeshType>().cellVolumes();

    Sys.A() = factor*cv/this->deltaT();
    Sys.b() = factor*cv*field.oldTime()/this->deltaT();

    if (lambdaPtr)
    {
        // The restricted lambda field has undefined ghosts. Manually implement
        // operation to avoid segmentation errors.

        forAllCells(Sys.A(), l, d, i, j, k)
            Sys.A()[l][d](i,j,k) *=
                (*lambdaPtr)[l][d](i,j,k);

        const meshField<scalar,MeshType>& lambdaOld = lambdaPtr->oldTime();

        forAllCells(Sys.b(), l, d, i, j, k)
            Sys.b()[l][d](i,j,k) *= lambdaOld[l][d](i,j,k);
    }

    if (lambdaPtr)
        const_cast<meshField<scalar,MeshType>&>(*lambdaPtr)
       .makeShallow();

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
EulerDdtScheme<Type,MeshType>::exDdt
(
    const meshField<scalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field,
    const scalar factor
)
{
    tmp<meshField<Type,MeshType>> tDdt
    (
        new meshField<Type,MeshType>
        (
            lambdaPtr
          ? "ddt("+lambdaPtr->name()+","+field.name()+")"
          : "ddt("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<Type,MeshType>& ddt = tDdt.ref();

    if (lambdaPtr)
    {
        ddt =
            factor
          * (
                (*lambdaPtr)*field
              - lambdaPtr->oldTime()*field.oldTime()
            )
          / this->deltaT();
    }
    else
    {
        ddt = factor*(field-field.oldTime())/this->deltaT();
    }

    return tDdt;
}

}

}

}
