#include "EulerDdtScheme.H"
#include "restrict.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
EulerDdtScheme<Type,MeshType>::EulerDdtScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    ddtScheme<Type,MeshType>(fvMsh, is)
{}

template<class Type, class MeshType>
tmp<linearSystem<diagStencil,Type,MeshType>>
EulerDdtScheme<Type,MeshType>::imDdt
(
    const meshField<scalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field
)
{
    if (lambdaPtr)
        restrict(*lambdaPtr);

    tmp<linearSystem<diagStencil,Type,MeshType>> tSys =
        linearSystem<diagStencil,Type,MeshType>::New
        (
            lambdaPtr
          ? word("ddt("+lambdaPtr->name()+","+field.name()+")")
          : word("ddt("+field.name()+")"),
            const_cast<meshField<Type,MeshType>&>(field)
        );

    linearSystem<diagStencil,Type,MeshType>& Sys = tSys.ref();

    const meshField<scalar,MeshType>& cv =
        this->fvMsh().template metrics<MeshType>().cellVolumes();

    Sys.A() = cv/this->deltaT();
    Sys.b() = cv*field.oldTime()/this->deltaT();

    if (lambdaPtr)
    {
        // The restricted lambda field has undefined ghosts. Manually implement
        // operation to avoid segmentation errors.

        forAllCells(Sys.A(), l, d, i, j, k)
            Sys.A()[l][d](i,j,k) *= (*lambdaPtr)[l][d](i,j,k);

        const meshField<scalar,MeshType>& lambdaOld = lambdaPtr->oldTime();

        forAllCells(Sys.b(), l, d, i, j, k)
            Sys.b()[l][d](i,j,k) *= lambdaOld[l][d](i,j,k);
    }

    if (lambdaPtr)
        collapse(*lambdaPtr);

    return tSys;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>>
EulerDdtScheme<Type,MeshType>::exDdt
(
    const meshField<scalar,MeshType>* lambdaPtr,
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<Type,MeshType>> tDdt =
        meshField<Type,MeshType>::New
        (
            lambdaPtr
          ? "ddt("+lambdaPtr->name()+","+field.name()+")"
          : "ddt("+field.name()+")",
            field.fvMsh()
        );

    meshField<Type,MeshType>& ddt = tDdt.ref();

    if (lambdaPtr)
    {
        ddt =
            (
                (*lambdaPtr)*field
              - lambdaPtr->oldTime()*field.oldTime()
            )
          / this->deltaT();
    }
    else
    {
        ddt = field-field.oldTime()/this->deltaT();
    }

    return tDdt;
}

}

}

}
