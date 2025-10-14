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
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
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
    const fvMesh& fvMsh,
    const word schemeName
)
{
    Istream& is = getStream(fvMsh, "divergenceSchemes", schemeName);

    word divergenceSchemeType;
    is >> divergenceSchemeType;

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(divergenceSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown divergence scheme "
            << divergenceSchemeType << nl << nl
            << "Valid divergence schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<divergenceScheme<SType,Type,MeshType>>
    (
        cstrIter()(fvMsh, is)
    );
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> div
(
    const meshField<FaceSpace<Type>,MeshType>& phi
)
{
    tmp<meshField<Type,MeshType>> tDiv =
        meshField<Type,MeshType>::New
        (
            "div("+phi.name()+")",
            phi.fvMsh()
        );

    meshField<Type,MeshType>& Div = tDiv.ref();

    const meshField<scalar,MeshType>& icv =
        phi.fvMsh().template metrics<MeshType>().inverseCellVolumes();

    forAllCells(phi, d, i, j, k)
    {
        #ifdef NO_BLOCK_ZERO_INIT
        Div(d,i,j,k) = Zero;
        #endif

        for (int f = 0; f < 6; f++)
            Div(d,i,j,k) += phi(d,i,j,k)[f];

        Div(d,i,j,k) *= icv(d,i,j,k);
    }

    return tDiv;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> div
(
    const tmp<meshField<FaceSpace<Type>,MeshType>>& tPhi
)
{
    tmp<meshField<Type,MeshType>> tDiv = div(tPhi());

    if (tPhi.isTmp())
        tPhi.clear();

    return tDiv;
}

template<class Type>
tmp<meshField<Type,colocated>> coloDiv
(
    const meshField<Type,staggered>& field
)
{
    tmp<meshField<Type,colocated>> tDiv =
        meshField<Type,colocated>::New
        (
            "coloDiv("+field.name()+")",
            field.fvMsh()
        );

    meshField<Type,colocated>& Div = tDiv.ref();

    const meshField<scalar,colocated>& icv =
        field.fvMsh().template metrics<colocated>().inverseCellVolumes();

    const meshField<faceScalar,colocated>& fa =
        field.fvMsh().template metrics<colocated>().faceAreas();

    forAllCells(Div, i, j, k)
    {
        #ifdef NO_BLOCK_ZERO_INIT
        Div(i,j,k) = Zero;
        #endif

        for (int fd = 0; fd < 3; fd++)
            Div(i,j,k) -=
                field(fd,i,j,k)*fa(i,j,k)[fd*2]
              - field(fd,upperNeighbor(i,j,k,fd))*fa(i,j,k)[fd*2+1];

        Div(i,j,k) *= icv(i,j,k);
    }

    return tDiv;
}

template<class Type>
tmp<meshField<Type,colocated>> coloDiv
(
    const tmp<meshField<Type,staggered>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<Type,colocated>> tColoDiv = coloDiv(tField());

    if (tField.isTmp())
        tField.clear();

    return tColoDiv;
}

}

}

}
