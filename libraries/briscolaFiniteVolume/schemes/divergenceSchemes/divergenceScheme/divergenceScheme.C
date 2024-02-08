#include "divergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
divergenceScheme<SType,Type,MeshType>::divergenceScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
{}

template<class SType, class Type, class MeshType>
divergenceScheme<SType,Type,MeshType>::divergenceScheme
(
    const divergenceScheme<SType,Type,MeshType>& s
)
:
    scheme(s)
{}

template<class SType, class Type, class MeshType>
divergenceScheme<SType,Type,MeshType>::~divergenceScheme()
{}

template<class SType, class Type, class MeshType>
autoPtr<divergenceScheme<SType,Type,MeshType>>
divergenceScheme<SType,Type,MeshType>::New
(
    const word name,
    const fvMesh& fvMsh
)
{
    const dictionary dict
    (
        fvMsh.schemeDict().subDict("divergenceSchemes").subDict(name)
    );

    const word divergenceSchemeType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(divergenceSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown divergence scheme "
            << divergenceSchemeType << nl << nl
            << "Valid divergence schemes are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<divergenceScheme<SType,Type,MeshType>>
    (
        cstrIter()(dict, fvMsh)
    );
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> explicitDiv
(
    const meshField<LowerFaceSpace<Type>,MeshType>& phi
)
{
    tmp<meshField<Type,MeshType>> tDiv
    (
        new meshField<Type,MeshType>
        (
            "div("+phi.name()+")",
            phi.fvMsh()
        )
    );

    meshField<Type,MeshType>& Div = tDiv.ref();

    Div = Zero;

    const meshField<scalar,MeshType>& cv =
        phi.fvMsh().template metrics<MeshType>().cellVolumes();

    forAllCells(phi, d, i, j, k)
        for (int fd = 0; fd < 3; fd++)
            Div(d,i,j,k) +=
                (
                    phi(d,i,j,k)[fd]
                  - phi(d,upperNei(i,j,k,fd))[fd]
                )
              / cv(d,i,j,k);

    return tDiv;
}

template<class Type>
tmp<meshField<Type,colocated>> explicitColoDiv
(
    const meshField<Type,staggered>& field
)
{
    tmp<meshField<Type,colocated>> tDiv
    (
        new meshField<Type,colocated>
        (
            "coloDiv("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<Type,colocated>& Div = tDiv.ref();

    Div = Zero;

    const meshField<scalar,colocated>& cv =
        field.fvMsh().template metrics<colocated>().cellVolumes();

    const meshField<faceScalar,colocated>& fa =
        field.fvMsh().template metrics<colocated>().faceAreas();

    forAllCells(Div, i, j, k)
        for (int fd = 0; fd < 3; fd++)
            Div(i,j,k) -=
                (
                    field(fd,i,j,k)*fa(i,j,k)[fd*2]
                  - field(fd,upperNei(i,j,k,fd))*fa(i,j,k)[fd*2+1]
                )
              / cv(i,j,k);

    return tDiv;
}

}

}

}
