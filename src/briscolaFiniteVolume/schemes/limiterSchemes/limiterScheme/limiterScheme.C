#include "limiterScheme.H"
#include "colocated.H"
#include "staggered.H"

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
    const meshField<scalar,MeshType>& field,
    const bool deep
)
{
    tmp<faceField<scalar,MeshType>> tR =
        faceField<scalar,MeshType>::New("r", field.fvMsh());

    faceField<scalar,MeshType>& R = tR.ref();

    R.make(deep);

    // Field

    if (sGradPtr_.empty())
        sGradPtr_.set
        (
            meshField<vector,MeshType>::New
            (
                "grad(" + field.name() + ")",
                field.fvMsh()
            ).ptr()
        );

    meshField<vector,MeshType>& grad = sGradPtr_();

    // Gradient

    if (sGradSchemePtr_.empty())
        sGradSchemePtr_.set
        (
            new midPointGaussGradientScheme<scalar,MeshType>
            (
                field.fvMsh()
            )
        );

    midPointGaussGradientScheme<scalar,MeshType>& gradScheme =
        sGradSchemePtr_();

    grad = gradScheme.grad(field);

    // Restriction

    if (deep)
    {
        if (sRestrictionSchemePtr_.empty())
            sRestrictionSchemePtr_.set
            (
                restrictionScheme<vector,MeshType>::New
                (
                    field.fvMsh(),
                    "average"
                ).ptr()
            );

        sRestrictionSchemePtr_->restrict(grad);
    }

    const meshField<vector,MeshType>& cc =
        field.fvMsh().template metrics<MeshType>().cellCenters();

    // Lambda function definition

    auto calc =
        [&](label fd, label l, label d, label i, label j, label k)
        -> scalar
    {
        // OpenFOAM's way (see NVDTVD.H)

        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const vector dist = cc(l,d,nei) - cc(l,d,ijk);
        const scalar gradf = field(l,d,nei) - field(l,d,ijk);

        scalar gradc;

        if (phi[fd](l,d,ijk) > 0)
        {
            gradc = dist & grad(l,d,ijk);
        }
        else
        {
            gradc = dist & grad(l,d,nei);
        }

        if (mag(gradc) >= 1000*mag(gradf))
        {
            return 2*1000*sign(gradc)*sign(gradf) - 1;
        }
        else
        {
            return 2*gradc/gradf - 1;
        }
    };

    // Communicate/compute

    const label nReq = Pstream::nRequests();

    grad.template prepare<bcsOfType<parallelBoundary>>();
    grad.correctUnsetBoundaryConditions();

    forAllInternalFaces(R, fd, l, d, i, j, k)
        R[fd](l,d,i,j,k) = calc(fd,l,d,i,j,k);

    if (Pstream::parRun())
        UPstream::waitRequests(nReq);

    grad.template evaluate<bcsOfType<parallelBoundary>>();

    forAllBoundaryFaces(R, fd, lu, l, d, i, j, k)
        R[fd](l,d,i,j,k) = calc(fd,l,d,i,j,k);

    return tR;
}

template<class Type, class MeshType>
tmp<faceField<scalar,MeshType>> limiterScheme<Type,MeshType>::rType
(
    const faceField<scalar,MeshType>& phi,
    const meshField<vector,MeshType>& field,
    const bool deep
)
{
    tmp<faceField<scalar,MeshType>> tR =
        faceField<scalar,MeshType>::New("r", field.fvMsh());

    faceField<scalar,MeshType>& R = tR.ref();

    R.make(deep);

    // Field

    if (vGradPtr_.empty())
        vGradPtr_.set
        (
            meshField<tensor,MeshType>::New
            (
                "grad(" + field.name() + ")",
                field.fvMsh()
            ).ptr()
        );

    meshField<tensor,MeshType>& grad = vGradPtr_();

    // Gradient

    if (vGradSchemePtr_.empty())
        vGradSchemePtr_.set
        (
            new midPointGaussGradientScheme<vector,MeshType>
            (
                field.fvMsh()
            )
        );

    grad = vGradSchemePtr_->grad(field);

    // Restriction

    if (deep)
    {
        if (vRestrictionSchemePtr_.empty())
            vRestrictionSchemePtr_.set
            (
                restrictionScheme<tensor,MeshType>::New
                (
                    field.fvMsh(),
                    "average"
                ).ptr()
            );

        vRestrictionSchemePtr_->restrict(grad);
    }

    const meshField<vector,MeshType>& cc =
        field.fvMsh().template metrics<MeshType>().cellCenters();

    // Lambda function definition

    auto calc =
        [&](label fd, label l, label d, label i, label j, label k)
        -> scalar
    {
        // OpenFOAM's way (see NVDTVDV.H)

        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const vector dist = cc(l,d,nei) - cc(l,d,ijk);

        const vector gradfv = field(l,d,nei) - field(l,d,ijk);
        const scalar gradf = gradfv & gradfv;

        scalar gradc;

        if (phi[fd](l,d,ijk) > 0)
        {
            gradc = gradfv & (dist & grad(l,d,ijk));
        }
        else
        {
            gradc = gradfv & (dist & grad(l,d,nei));
        }

        if (mag(gradc) >= 1000*mag(gradf))
        {
           return 2*1000*sign(gradc)*sign(gradf) - 1;
        }
        else
        {
            return 2*gradc/gradf - 1;
        }
    };

    // Communicate/compute

    const label nReq = Pstream::nRequests();

    grad.template prepare<bcsOfType<parallelBoundary>>();
    grad.correctUnsetBoundaryConditions();

    forAllInternalFaces(R, fd, l, d, i, j, k)
        R[fd](l,d,i,j,k) = calc(fd,l,d,i,j,k);

    if (Pstream::parRun())
        UPstream::waitRequests(nReq);

    grad.template evaluate<bcsOfType<parallelBoundary>>();

    forAllBoundaryFaces(R, fd, lu, l, d, i, j, k)
        R[fd](l,d,i,j,k) = calc(fd,l,d,i,j,k);

    return tR;
}

}

}

}
