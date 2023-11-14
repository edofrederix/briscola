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
    interpolationScheme<Type,MeshType>(dict,fvMsh)
{}

template<class Type, class MeshType>
curvatureInterpolationScheme<Type,MeshType>::curvatureInterpolationScheme
(
    const fvMesh& fvMsh
)
:
    interpolationScheme<Type,MeshType>(dictionary(),fvMsh)
{}

template<class Type, class MeshType>
tmp<meshField<FaceSpace<Type>,MeshType>>
curvatureInterpolationScheme<Type,MeshType>::interp
(
    const meshField<Type,MeshType>& field
)
{
    tmp<meshField<FaceSpace<Type>,MeshType>> tInterp
    (
        new meshField<FaceSpace<Type>,MeshType>
        (
            "interp("+field.name()+")",
            field.fvMsh()
        )
    );

    meshField<FaceSpace<Type>,MeshType>& Interp = tInterp.ref();

    Interp = Zero;

    colocatedScalarField alpha = field.fvMsh().db().template
            lookupObject<colocatedScalarField>("alpha");

    forAllCells(field, i, j, k)
    {
        const labelVector ijk(i,j,k);
        Interp(ijk) = Zero;

        if
        (
            alpha(ijk) > 1e-12
            && alpha(ijk) < (1.0 - 1e-12)
        )
        {
            for (int dir = 0; dir < 2; dir++)
            {
                const labelVector ijku(ijk + units[dir]);
                const labelVector ijkl(ijk - units[dir]);

                if
                (
                    alpha(ijku) > 1e-12
                    && alpha(ijku) < (1.0 - 1e-12)
                )
                {
                    Interp(ijk)[2*dir+1] = 0.5 * (field(ijk) + field(ijku));
                }
                else
                {
                    Interp(ijk)[2*dir+1] = field(ijk);
                }

                if
                (
                    alpha(ijkl) > 1e-12
                    && alpha(ijkl) < (1.0 - 1e-12)
                )
                {
                    Interp(ijk)[2*dir] = 0.5 * (field(ijk) + field(ijkl));
                }
                else
                {
                    Interp(ijk)[2*dir] = field(ijk);
                }

            }

        }
        else
        {
            for (int dir = 0; dir < 2; dir++)
            {
                const labelVector ijku(ijk + units[dir]);
                const labelVector ijkl(ijk - units[dir]);

                Interp(ijk)[2*dir+1] = field(ijku);
                Interp(ijk)[2*dir] = field(ijkl);

            }
        }
    }

    return tInterp;
}

}

}

}
