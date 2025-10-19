#include "limiterScheme.H"
#include "colocated.H"
#include "staggered.H"
#include "exSchemesGradient.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
limiterScheme<Type,MeshType>::limiterScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    scheme(fvMsh)
{}

template<class Type, class MeshType>
limiterScheme<Type,MeshType>::limiterScheme
(
    const limiterScheme<Type,MeshType>& s
)
:
    scheme(s)
{}

template<class Type, class MeshType>
limiterScheme<Type,MeshType>::~limiterScheme()
{}

template<class Type, class MeshType>
autoPtr<limiterScheme<Type,MeshType>>
limiterScheme<Type,MeshType>::New
(
    const fvMesh& fvMsh,
    const word limiterSchemeType
)
{
    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(limiterSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown limiter scheme "
            << limiterSchemeType << nl << nl
            << "Valid limiter schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<limiterScheme<Type,MeshType>>
    (
        cstrIter()(fvMsh, scheme::nullStream)
    );
}

template<class Type, class MeshType>
autoPtr<limiterScheme<Type,MeshType>>
limiterScheme<Type,MeshType>::New
(
    const fvMesh& fvMsh,
    Istream& is
)
{
    word limiterSchemeType;
    is >> limiterSchemeType;

    typename IstreamConstructorTable::iterator cstrIter =
        IstreamConstructorTablePtr_->find(limiterSchemeType);

    if (cstrIter == IstreamConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown limiter scheme "
            << limiterSchemeType << nl << nl
            << "Valid limiter schemes are:" << nl
            << IstreamConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<limiterScheme<Type,MeshType>>
    (
        cstrIter()(fvMsh, is)
    );
}

template<class Type, class MeshType>
tmp<faceField<scalar,MeshType>> limiterScheme<Type,MeshType>::rType
(
    const faceField<scalar,MeshType>& phi,
    const meshField<scalar,MeshType>& field
)
{
    tmp<faceField<scalar,MeshType>> tR =
        faceField<scalar,MeshType>::New("r", field.fvMsh());

    faceField<scalar,MeshType>& R = tR.ref();

    forAll(R, fd)
        R[fd].make(phi[fd].deep() && field.deep());

    meshField<vector,MeshType> grad
    (
        gradientScheme<scalar,MeshType>::New
        (
            field.fvMsh(),
            "grad("+field.name()+")"
        )->grad(field)
    );

    grad.correctBoundaryConditions();

    if (R[0].deep())
        restrict(grad);

    const meshField<vector,MeshType>& cc =
        field.fvMsh().template metrics<MeshType>().cellCenters();

    forAllFaces(R, fd, l, d, i, j, k)
    {
        // OpenFOAM's way (see NVDTVD.H)

        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const vector dist = cc(l,d,nei) - cc(l,d,ijk);
        const scalar gradf = field(l,d,nei) - field(l,d,ijk);
        const scalar gradc =
            dist & (phi[fd](l,d,ijk) > 0 ? grad(l,d,ijk) : grad(l,d,nei));

        R[fd](l,d,ijk) =
            Foam::mag(gradc) >= 1000*Foam::mag(gradf)
          ? 2*1000*Foam::sign(gradc)*Foam::sign(gradf) - 1
          : 2*(gradc/gradf) - 1;
    }

    return tR;
}

template<class Type, class MeshType>
tmp<faceField<scalar,MeshType>> limiterScheme<Type,MeshType>::rType
(
    const faceField<scalar,MeshType>& phi,
    const meshField<vector,MeshType>& field
)
{
    tmp<faceField<scalar,MeshType>> tR =
        faceField<scalar,MeshType>::New("r", field.fvMsh());

    faceField<scalar,MeshType>& R = tR.ref();

    forAll(R, fd)
        R[fd].make(phi[fd].deep() && field.deep());

    meshField<tensor,MeshType> grad
    (
        gradientScheme<vector,MeshType>::New
        (
            field.fvMsh(),
            "grad("+field.name()+")"
        )->grad(field)
    );

    grad.correctBoundaryConditions();

    if (R[0].deep())
        restrict(grad);

    const meshField<vector,MeshType>& cc =
        field.fvMsh().template metrics<MeshType>().cellCenters();

    forAllFaces(R, fd, l, d, i, j, k)
    {
        // OpenFOAM's way (see NVDTVDV.H)

        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const vector dist = cc(l,d,nei) - cc(l,d,ijk);

        const vector gradfv = field(l,d,nei) - field(l,d,ijk);
        const scalar gradf = gradfv & gradfv;

        const scalar gradc =
            gradfv
          & (
                dist
              & (
                    phi[fd](l,d,ijk) > 0
                  ? grad(l,d,ijk)
                  : grad(l,d,nei)
                )
            );

        R[fd](l,d,ijk) =
            Foam::mag(gradc) >= 1000*Foam::mag(gradf)
          ? 2*1000*Foam::sign(gradc)*Foam::sign(gradf) - 1
          : 2*(gradc/gradf) - 1;
    }

    return tR;
}

}

}

}
