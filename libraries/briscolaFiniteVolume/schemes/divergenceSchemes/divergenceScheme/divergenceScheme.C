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

    forAll(phi, l)
    forAll(phi[l], d)
    {
        meshDirection<Type,MeshType>& D = Div[l][d];
        const meshDirection<FaceSpace<Type>,MeshType>& p = phi[l][d];

        const meshDirection<scalar,MeshType>& cv =
            phi.fvMsh().template
            metrics<MeshType>().cellVolumes()[l][d];

        D.initGhosts();

        forAllCells(p, i, j, k)
        {
            D(i,j,k) = neighborSum(p(i,j,k))/cv(i,j,k);
        }
    }

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

    forAll(Div, l)
    {
        const meshDirection<Type,staggered>& f0 = field[l][0];
        const meshDirection<Type,staggered>& f1 = field[l][1];
        const meshDirection<Type,staggered>& f2 = field[l][2];

        meshDirection<Type,colocated>& D = Div[l][0];

        const meshDirection<scalar,colocated>& cv =
            field.fvMsh().template
            metrics<colocated>().cellVolumes()[l][0];

        const meshDirection<faceScalar,colocated>& fa =
            field.fvMsh().template
            metrics<colocated>().faceAreas()[l][0];

        D.initGhosts();

        forAllCells(D, i, j, k)
        {
            D(i,j,k) =
                (
                  - f0(i,  j,  k  ) * fa(i,j,k).left()
                  + f0(i+1,j,  k  ) * fa(i,j,k).right()
                  - f1(i,  j,  k  ) * fa(i,j,k).bottom()
                  + f1(i,  j+1,k  ) * fa(i,j,k).top()
                  - f2(i,  j,  k  ) * fa(i,j,k).aft()
                  + f2(i,  j,  k+1) * fa(i,j,k).fore()
                )
              / cv(i,j,k);
        }
    }

    return tDiv;
}

}

}

}
