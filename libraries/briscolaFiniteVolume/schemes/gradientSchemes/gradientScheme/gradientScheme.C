#include "gradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
gradientScheme<Type,MeshType>::gradientScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    scheme(dict, fvMsh)
{}

template<class Type, class MeshType>
gradientScheme<Type,MeshType>::gradientScheme
(
    const gradientScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
gradientScheme<Type,MeshType>::~gradientScheme()
{}

template<class Type, class MeshType>
autoPtr<gradientScheme<Type,MeshType>>
gradientScheme<Type,MeshType>::New
(
    const word name,
    const fvMesh& fvMsh
)
{
    const dictionary dict
    (
        fvMsh.schemeDict().subDict("gradientSchemes").subDict(name)
    );

    const word gradientSchemeType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(gradientSchemeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown gradient scheme "
            << gradientSchemeType << nl << nl
            << "Valid gradient schemes are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<gradientScheme<Type,MeshType>>
    (
        cstrIter()(dict, fvMsh)
    );
}

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
centerGrad
(
    const meshField<FaceSpace<Type>,MeshType>& field
)
{
    tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>> tGrad
    (
        new meshField<typename outerProduct<vector,Type>::type,MeshType>
        (
            "centerGrad("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<typename outerProduct<vector,Type>::type,MeshType>& Grad =
        tGrad.ref();

    Grad = Zero;

    const meshField<scalar,MeshType>& cv =
        field.fvMsh().template metrics<MeshType>().cellVolumes();

    const meshField<faceVector,MeshType>& fan =
        field.fvMsh().template metrics<MeshType>().faceAreaNormals();

    forAllCells(Grad, d, i, j, k)
        for (int f = 0; f < 6; f++)
            Grad(d,i,j,k) +=
                field(d,i,j,k)[f]*fan(d,i,j,k)[f]/cv(d,i,j,k);

    return tGrad;
}

}

}

}
