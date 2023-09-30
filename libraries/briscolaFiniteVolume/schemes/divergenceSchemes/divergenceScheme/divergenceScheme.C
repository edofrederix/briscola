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

    forAllLevels(Div, l, d, i, j, k)
        Div(l,d,i,j,k) = neighborSum(phi(l,d,i,j,k))/cv(l,d,i,j,k);

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

    forAllLevels(Div, l, d, i, j, k)
        Div(l,d,i,j,k) =
            (
              - field(l,0,i,  j,  k  ) * fa(l,d,i,j,k).left()
              + field(l,0,i+1,j,  k  ) * fa(l,d,i,j,k).right()
              - field(l,1,i,  j,  k  ) * fa(l,d,i,j,k).bottom()
              + field(l,1,i,  j+1,k  ) * fa(l,d,i,j,k).top()
              - field(l,2,i,  j,  k  ) * fa(l,d,i,j,k).aft()
              + field(l,2,i,  j,  k+1) * fa(l,d,i,j,k).fore()
            )
          / cv(l,d,i,j,k);

    return tDiv;
}

}

}

}
