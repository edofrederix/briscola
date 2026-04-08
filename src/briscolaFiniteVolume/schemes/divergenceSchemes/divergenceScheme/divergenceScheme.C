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
    const faceField<Type,MeshType>& phi
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

    Div = Zero;

    forAllFaces(phi, fd, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        Div(d,ijk) += phi[fd](d,ijk);
        Div(d,nei) -= phi[fd](d,ijk);
    }

    Div *= icv;

    return tDiv;
}

template<class Type, class MeshType>
tmp<meshField<Type,MeshType>> div
(
    const tmp<faceField<Type,MeshType>>& tPhi
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

    const faceField<scalar,colocated>& fa =
        field.fvMsh().template metrics<colocated>().faceAreas();

    Div = Zero;

    forAllFaces(fa, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const Type value = field(fd,ijk)*fa[fd](ijk);

        Div(ijk) -= value;
        Div(nei) += value;
    }

    Div *= icv;

    return tDiv;
}

template<class Type>
tmp<meshField<Type,colocated>> coloDiv
(
    const tmp<meshField<Type,staggered>>& tField
)
{
    if (tField.isTmp())
        const_cast<tmp<meshField<Type,staggered>>&>(tField)
            ->correctBoundaryConditions();

    tmp<meshField<Type,colocated>> tColoDiv = coloDiv(tField());

    if (tField.isTmp())
        tField.clear();

    return tColoDiv;
}

}

}

}
