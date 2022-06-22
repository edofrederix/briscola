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
EulerDdtScheme<Type,MeshType>::ddt
(
    meshField<Type,MeshType>& field
)
{
    tmp<linearSystem<diagStencil,Type,MeshType>> tSys
    (
        new linearSystem<diagStencil,Type,MeshType>(field)
    );

    linearSystem<diagStencil,Type,MeshType>& Sys = tSys.ref();

    const meshField<scalar,MeshType>& cv =
        this->fvMsh().template metrics<MeshType>().cellVolumes();

    Sys.A() = cv/this->deltaT();
    Sys.b() = cv*field.oldTime()/this->deltaT();

    return tSys;
}

template<class Type, class MeshType>
tmp<linearSystem<diagStencil,Type,MeshType>>
EulerDdtScheme<Type,MeshType>::ddt
(
    const meshField<scalar,MeshType>& coeff,
    meshField<Type,MeshType>& field
)
{
    tmp<linearSystem<diagStencil,Type,MeshType>> tSys
    (
        new linearSystem<diagStencil,Type,MeshType>(field)
    );

    linearSystem<diagStencil,Type,MeshType>& Sys = tSys.ref();

    const meshField<scalar,MeshType>& cv =
        this->fvMsh().template metrics<MeshType>().cellVolumes();

    Sys.A() = cv*coeff/this->deltaT();
    Sys.b() = cv*coeff.oldTime()*field.oldTime()/this->deltaT();

    return tSys;
}

}

}

}
