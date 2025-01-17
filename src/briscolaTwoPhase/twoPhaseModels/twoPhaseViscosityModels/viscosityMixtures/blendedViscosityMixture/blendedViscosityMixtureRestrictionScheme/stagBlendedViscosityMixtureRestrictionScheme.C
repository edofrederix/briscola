#include "stagBlendedViscosityMixtureRestrictionScheme.H"

#include "restrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"
#include "staggered.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeRestrictionSchemeNoTemplate
(
    stagBlendedViscosityMixture,
    faceScalar,
    staggered
);

defineTemplateTypeNameAndDebugWithName
(
    blendedViscosityMixtureRestrictionScheme<staggered>,
    "blendedViscosityMixture",
    0
);

stagBlendedViscosityMixtureRestrictionScheme::
stagBlendedViscosityMixtureRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    blendedViscosityMixtureRestrictionScheme<staggered>(fvMsh, is)
{}

void stagBlendedViscosityMixtureRestrictionScheme::restrict
(
    meshDirection<faceScalar,staggered>& coarse,
    const meshDirection<faceScalar,staggered>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());
    const label d = coarse.directionNum();

    forAllFaces(coarse, fd, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const labelVector nei(lowerNeighbor(i,j,k,fd));

        const labelVector fijk(briscola::cmptMultiply(ijk,R));

        labelVector o(Zero);

        if (fd == d)
            o[d] = -1;

        // Reconstruct the face alpha

        scalar alpha = 0.0;
        label nFaces = 0;

        for (o.x(); o.x() < (fd == 0 || d == 0 ? 1 : R.x()); o.x()++)
        for (o.y(); o.y() < (fd == 1 || d == 1 ? 1 : R.y()); o.y()++)
        for (o.z(); o.z() < (fd == 2 || d == 2 ? 1 : R.z()); o.z()++)
        {
            // Don't use cells with negative indices
            const labelVector fijko(briscola::cmptMax(fijk+o,zeroXYZ));

            alpha += this->inv(fine(fijko)[fd*2]);
            nFaces++;
        }

        alpha /= scalar(nFaces);

        // Apply the blending function. Viscosity has no sign so we can copy
        // directly to the upper face of the lower neighbor.

        coarse(ijk)[fd*2  ] = this->blend(alpha);
        coarse(ijk)[fd*2+1] = coarse(ijk)[fd*2];
    }
}

}

}

}
