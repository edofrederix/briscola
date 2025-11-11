#include "limiterScheme.H"
#include "colocated.H"
#include "staggered.H"
#include "midPointGaussGradientScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

using ::Foam::mag;
using ::Foam::sign;

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

    tmp<meshField<vector,MeshType>> tGrad =
        midPointGaussGradientScheme<scalar,MeshType>(field.fvMsh()).grad(field);

    meshField<vector,MeshType>& grad = tGrad.ref();
    grad.correctCommsBoundaryConditions();

    const meshField<vector,MeshType>& cc =
        field.fvMsh().template metrics<MeshType>().cellCenters();

    forAllFaces(R, fd, d, i, j, k)
    {
        // OpenFOAM's way (see NVDTVD.H)

        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const vector dist = cc(d,nei) - cc(d,ijk);
        const scalar gradf = field(d,nei) - field(d,ijk);

        scalar gradc;

        if (phi[fd](d,ijk) > 0)
        {
            gradc = dist & grad(d,ijk);
        }
        else
        {
            gradc = dist & grad(d,nei);
        }

        if (mag(gradc) >= 1000*mag(gradf))
        {
            R[fd](d,ijk) = 2*1000*sign(gradc)*sign(gradf) - 1;
        }
        else
        {
            R[fd](d,ijk) = 2*gradc/gradf - 1;
        }
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

    tmp<meshField<tensor,MeshType>> tGrad =
        midPointGaussGradientScheme<vector,MeshType>(field.fvMsh()).grad(field);

    meshField<tensor,MeshType>& grad = tGrad.ref();
    grad.correctCommsBoundaryConditions();

    const meshField<vector,MeshType>& cc =
        field.fvMsh().template metrics<MeshType>().cellCenters();

    forAllFaces(R, fd, d, i, j, k)
    {
        // OpenFOAM's way (see NVDTVDV.H)

        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const vector dist = cc(d,nei) - cc(d,ijk);

        const vector gradfv = field(d,nei) - field(d,ijk);
        const scalar gradf = gradfv & gradfv;

        scalar gradc;

        if (phi[fd](d,ijk) > 0)
        {
            gradc = gradfv & (dist & grad(d,ijk));
        }
        else
        {
            gradc = gradfv & (dist & grad(d,nei));
        }

        if (mag(gradc) >= 1000*mag(gradf))
        {
            R[fd](d,ijk) = 2*1000*sign(gradc)*sign(gradf) - 1;
        }
        else
        {
            R[fd](d,ijk) = 2*gradc/gradf - 1;
        }
    }

    return tR;
}

}

}

}
