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
    restrictionScheme<faceScalar,colocated>(fvMsh, is)
{}

void coloFluxRestrictionScheme::restrict
(
    meshDirection<faceScalar,colocated>& coarse,
    const meshDirection<faceScalar,colocated>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());

    // Sum fluxes from corresponding fine grid faces

    forAllFaces(coarse, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const labelVector fijk(briscola::cmptMultiply(ijk,R));
        const labelVector fnei(lowerNeighbor(fijk,fd));

        coarse(ijk)[fd*2  ] = Zero;
        coarse(nei)[fd*2+1] = Zero;

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            coarse(ijk)[fd*2  ] += fine(fijk+o)[fd*2  ];
            coarse(nei)[fd*2+1] += fine(fnei+o)[fd*2+1];
        }
    }
}

}

}

}
