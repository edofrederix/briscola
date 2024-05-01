#include "curvatureInterpolationScheme.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
curvatureInterpolationScheme<Type,MeshType>::curvatureInterpolationScheme
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    interpolationScheme<Type,MeshType>(dict,fvMsh),
    threshold_(dict.lookupOrDefault<scalar>("threshold", 1e-12))
{}

template<class Type, class MeshType>
curvatureInterpolationScheme<Type,MeshType>::curvatureInterpolationScheme
(
    const fvMesh& fvMsh
)
:
    interpolationScheme<Type,MeshType>(dictionary(),fvMsh),
    threshold_(dictionary().lookupOrDefault<scalar>("threshold", 1e-12))
{}

template<>
tmp<colocatedLowerFaceScalarField>
curvatureInterpolationScheme<scalar,colocated>::interp
(
    const colocatedScalarField& field
)
{
    tmp<colocatedLowerFaceScalarField> tInterp
    (
        new colocatedLowerFaceScalarField
        (
            "interp("+field.name()+")",
            field.fvMsh()
        )
    );

    colocatedLowerFaceScalarField& Interp = tInterp.ref();

    Interp = Zero;

    const colocatedScalarField& alpha = field.fvMsh().db().template
            lookupObjectRef<colocatedScalarField>("alpha");

    forAllCells(field, i, j, k)
    {
        const labelVector ijk(i,j,k);
        if
        (
            alpha(ijk) > threshold_
            && alpha(ijk) < (1.0 - threshold_)
        )
        {
            for (int dir = 0; dir < 3; dir++)
            {
                const labelVector ijku(ijk + units[dir]);
                const labelVector ijkl(ijk - units[dir]);

                if
                (
                    alpha(ijku) > threshold_
                 && alpha(ijku) < (1.0 - threshold_)
                )
                {
                    Interp(ijku)[dir] = 0.5 * (field(ijk) + field(ijku));
                }
                else
                {
                    Interp(ijku)[dir] = field(ijk);
                }

                if
                (
                    alpha(ijkl) > threshold_
                    && alpha(ijkl) < (1.0 - threshold_)
                )
                {
                    Interp(ijk)[dir] = 0.5 * (field(ijk) + field(ijkl));
                }
                else
                {
                    Interp(ijk)[dir] = field(ijk);
                }

            }

        }
        else
        {
            for (int dir = 0; dir < 3; dir++)
            {
                const labelVector ijku(ijk + units[dir]);
                const labelVector ijkl(ijk - units[dir]);

                Interp(ijku)[dir] = field(ijku);
                Interp(ijk)[dir] = field(ijkl);

            }
        }
    }

    return tInterp;
}

template<>
tmp<staggeredLowerFaceScalarField>
curvatureInterpolationScheme<scalar,staggered>::interp
(
    const staggeredScalarField& field
)
{
    tmp<staggeredLowerFaceScalarField> tInterp
    (
        new staggeredLowerFaceScalarField
        (
            "interp("+field.name()+")",
            field.fvMsh()
        )
    );

    staggeredLowerFaceScalarField& Interp = tInterp.ref();

    Interp = Zero;

    const colocatedScalarField& alphaColocated = field.fvMsh().db().template
            lookupObject<colocatedScalarField>("alpha");

    staggeredScalarField alpha
        (
            "alpha",
            field.fvMsh()
        );

    forAllCells(alpha, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector ijkl(ijk - units[d]);
        alpha(d,i,j,k) = 0.5 * (alphaColocated(ijk) + alphaColocated(ijkl));
    }

    forAllCells(field, d, i, j, k)
    {
        const labelVector ijk(i,j,k);
        Interp(ijk) = Zero;

        if
        (
            alpha(d, ijk) > threshold_
            && alpha(d, ijk) < (1.0 - threshold_)
        )
        {
            for (int dir = 0; dir < 3; dir++)
            {
                const labelVector ijku(ijk + units[dir]);
                const labelVector ijkl(ijk - units[dir]);

                if
                (
                    alpha(d, ijku) > threshold_
                    && alpha(d, ijku) < (1.0 - threshold_)
                )
                {
                    Interp(d, ijku)[dir] = 0.5 * (field(d, ijk) + field(d, ijku));
                }
                else
                {
                    Interp(d, ijku)[dir] = field(d, ijk);
                }

                if
                (
                    alpha(d, ijkl) > threshold_
                    && alpha(d, ijkl) < (1.0 - threshold_)
                )
                {
                    Interp(d, ijk)[dir] = 0.5 * (field(d, ijk) + field(d, ijkl));
                }
                else
                {
                    Interp(d, ijk)[dir] = field(d, ijk);
                }

            }

        }
        else
        {
            for (int dir = 0; dir < 3; dir++)
            {
                const labelVector ijku(ijk + units[dir]);
                const labelVector ijkl(ijk - units[dir]);

                Interp(d, ijku)[dir] = field(d, ijku);
                Interp(d, ijk)[dir] = field(d, ijkl);

            }
        }
    }

    return tInterp;
}

}

}

}
