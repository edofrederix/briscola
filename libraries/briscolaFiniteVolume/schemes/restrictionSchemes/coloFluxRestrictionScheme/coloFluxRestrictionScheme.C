#include "coloFluxRestrictionScheme.H"

#include "colocated.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

coloFluxRestrictionScheme::coloFluxRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    restrictionScheme<lowerFaceScalar,colocated>(fvMsh, is)
{}

void coloFluxRestrictionScheme::restrict
(
    meshDirection<lowerFaceScalar,colocated>& coarse,
    const meshDirection<lowerFaceScalar,colocated>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());

    coarse = Zero;

    // Sum fluxes from corresponding fine grid faces

    forAllFaces(coarse, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector nei(ijk-units[fd]);

        labelVector fijk(i*R.x(), j*R.y(), k*R.z());
        labelVector fnei(fijk-units[fd]);

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            coarse(ijk)[fd] += fine(fijk+o)[fd];
        }
    }
}

}

}

}
