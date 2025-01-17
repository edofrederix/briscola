#include "coloBlendedViscosityMixtureRestrictionScheme.H"

#include "restrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"
#include "colocated.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeRestrictionSchemeNoTemplate
(
    coloBlendedViscosityMixture,
    faceScalar,
    colocated
);

defineTemplateTypeNameAndDebugWithName
(
    blendedViscosityMixtureRestrictionScheme<colocated>,
    "blendedViscosityMixture",
    0
);

coloBlendedViscosityMixtureRestrictionScheme::
coloBlendedViscosityMixtureRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    blendedViscosityMixtureRestrictionScheme<colocated>(fvMsh, is)
{}

void coloBlendedViscosityMixtureRestrictionScheme::restrict
(
    meshDirection<faceScalar,colocated>& coarse,
    const meshDirection<faceScalar,colocated>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());

    const meshDirection<faceScalar,colocated>& faf =
        this->fvMsh().template
        metrics<colocated>().faceAreas()
        [fine.levelNum()][fine.directionNum()];

    forAllFaces(coarse, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const labelVector fijk(briscola::cmptMultiply(ijk,R));

        // Reconstruct the face alpha

        scalar alpha = 0.0;
        scalar area = 0.0;

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            alpha += faf(fijk+o)[fd*2]*this->inv(fine(fijk+o)[fd*2]);
            area +=  faf(fijk+o)[fd*2];
        }

        alpha /= area;

        // Apply the blending function. Viscosity has no sign so we can copy
        // directly to the upper face of the lower neighbor.

        coarse(ijk)[fd*2  ] = this->blend(alpha);
        coarse(nei)[fd*2+1] = coarse(ijk)[fd*2];
    }
}

}

}

}
