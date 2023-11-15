#include "divergenceScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
divergenceScheme<Type,MeshType>::divergenceScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
{}

template<class Type, class MeshType>
divergenceScheme<Type,MeshType>::divergenceScheme
(
    const divergenceScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
divergenceScheme<Type,MeshType>::~divergenceScheme()
{}

template<class Type, class MeshType>
autoPtr<divergenceScheme<Type,MeshType>> divergenceScheme<Type,MeshType>::New
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

    return autoPtr<divergenceScheme<Type,MeshType>>(cstrIter()(dict, fvMsh));
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> explicitDiv
(
    const meshField<FaceSpace<Type>,MeshType>& phi
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

    forAllDirections(Div, d, i, j, k)
        Div(d,i,j,k) = neighborSum(phi(d,i,j,k))/cv(d,i,j,k);

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
        Div(i,j,k) =
            (
              - field(0,i,  j,  k  ) * fa(i,j,k).left()
              + field(0,i+1,j,  k  ) * fa(i,j,k).right()
              - field(1,i,  j,  k  ) * fa(i,j,k).bottom()
              + field(1,i,  j+1,k  ) * fa(i,j,k).top()
              - field(2,i,  j,  k  ) * fa(i,j,k).aft()
              + field(2,i,  j,  k+1) * fa(i,j,k).fore()
            )
          / cv(i,j,k);

    return tDiv;
}

}

}

}
