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

    Grad = Zero;

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

template<class Type, class MeshType>
tmp<meshField<Type,staggered>>
stagGrad
(
    const meshField<Type,colocated>& field
)
{
    tmp<meshField<Type,staggered>> tGrad =
        meshField<Type,staggered>::New
        (
            "stagGrad("+field.name()+")",
            field.fvMsh()
        );

    meshField<Type,staggered>& grad = tGrad.ref();

    const faceField<scalar,colocated>& delta =
            field.fvMsh().template metrics<colocated>().faceDeltas();

    forAllCells(grad, d, i, j, k)
        grad(d,i,j,k) =
            (field(i,j,k) - field(lowerNeighbor(i,j,k,d)))*delta[d](i,j,k);

    return tGrad;
}

template<class Type>
tmp<meshField<Type,staggered>>
stagGrad
(
    const tmp<meshField<Type,colocated>>& tField
)
{
    if (tField.isTmp())
        tField->correctBoundaryConditions();

    tmp<meshField<Type,staggered>> tStagGrad = stagGrad(tField());

    if (tField.isTmp())
        tField.clear();

    return tStagGrad;
}

}

}

}
