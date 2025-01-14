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
        const labelVector fnei(briscola::cmptMultiply(nei,R));

        labelVector o(Zero);

        if (fd == d)
            o[d] = -1;

        // Reconstruct the face alpha

        scalar alphaLower = 0.0;
        scalar alphaUpper = 0.0;

        label nFaces = 0;

        for (o.x(); o.x() < (fd == 0 || d == 0 ? 1 : R.x()); o.x()++)
        for (o.y(); o.y() < (fd == 1 || d == 1 ? 1 : R.y()); o.y()++)
        for (o.z(); o.z() < (fd == 2 || d == 2 ? 1 : R.z()); o.z()++)
        {
            // Don't use cells with negative indices
            const labelVector fijko(briscola::cmptMax(fijk+o,zeroXYZ));
            const labelVector fneio(briscola::cmptMax(fnei+o,zeroXYZ));

            alphaLower += this->inv(fine(fijko)[fd*2  ]);
            alphaUpper += this->inv(fine(fneio)[fd*2+1]);

            nFaces++;
        }

        alphaLower /= scalar(nFaces);
        alphaUpper /= scalar(nFaces);

        // Apply the blending function

        coarse(ijk)[fd*2  ] = this->blend(alphaLower);
        coarse(ijk)[fd*2+1] = this->blend(alphaUpper);
    }
}

}

}

}
