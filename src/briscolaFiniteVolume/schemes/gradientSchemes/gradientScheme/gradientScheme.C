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
    const faceField<Type,MeshType>& field
)
{
    typedef typename outerProduct<vector,Type>::type TypeR;

    tmp<meshField<TypeR,MeshType>> tGrad =
        meshField<TypeR,MeshType>::New
        (
            "grad("+field.name()+")",
            field.fvMsh()
        );

    meshField<TypeR,MeshType>& Grad =
        tGrad.ref();

    const meshField<scalar,MeshType>& icv =
        field.fvMsh().template metrics<MeshType>().inverseCellVolumes();

    const faceField<vector,MeshType>& fan =
        field.fvMsh().template metrics<MeshType>().faceAreaNormals();

    #ifdef NO_BLOCK_ZERO_INIT
    Grad = Zero;
    #endif

    forAllFaces(fan, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const TypeR value =
            field[fd](d,ijk)*fan[fd](d,ijk);

        Grad(d,ijk) += value;
        Grad(d,nei) -= value;
    }

    Grad *= icv;

    return tGrad;
}

template<class Type, class MeshType>
tmp<meshField<typename outerProduct<vector,Type>::type,MeshType>>
grad
(
    const tmp<faceField<Type,MeshType>>& tField
)
{
    typedef typename outerProduct<vector,Type>::type TypeR;

    tmp<meshField<TypeR,MeshType>> tGrad = grad(tField());

    if (tField.isTmp())
        tField.clear();

    return tGrad;
}

}

}

}
