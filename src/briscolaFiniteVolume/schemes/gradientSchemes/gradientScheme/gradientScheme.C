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
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
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
    const fvMesh& fvMsh,
    const word schemeName
)
{
    Istream& is = getStream(fvMsh, "gradientSchemes", schemeName);

    word gradientSchemeType;
    is >> gradientSchemeType;

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(gradientSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown gradient scheme "
            << gradientSchemeType << nl << nl
            << "Valid gradient schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<gradientScheme<Type,MeshType>>
    (
        cstrIter()(fvMsh, is)
    );
}

template<class Type, class MeshType>
inline tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const meshField<FaceSpace<Type>,MeshType>& field
)
{
    tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>> tGrad
    (
        new meshField<typename outerProduct<vector,Type>::type,MeshType>
        (
            "grad("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<typename outerProduct<vector,Type>::type,MeshType>& Grad =
        tGrad.ref();

    const meshField<scalar,MeshType>& cv =
        field.fvMsh().template metrics<MeshType>().cellVolumes();

    const meshField<faceVector,MeshType>& fan =
        field.fvMsh().template metrics<MeshType>().faceAreaNormals();

    forAllCells(Grad, d, i, j, k)
    {
        #ifdef NO_BLOCK_ZERO_INIT
        Grad(d,i,j,k) = Zero;
        #endif

        for (int fd = 0; fd < 3; fd++)
            Grad(d,i,j,k) +=
                (
                    field(d,i,j,k)[fd*2  ]*fan(d,i,j,k)[fd*2  ]
                  + field(d,i,j,k)[fd*2+1]*fan(d,i,j,k)[fd*2+1]
                )
              / cv(d,i,j,k);
    }

    return tGrad;
}

template<class Type, class MeshType>
tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const tmp<meshField<FaceSpace<Type>,MeshType>>& tField
)
{
    tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>> tGrad
    (
        grad(tField())
    );

    if (tField.isTmp())
        tField.clear();

    return tGrad;
}

}

}

}
