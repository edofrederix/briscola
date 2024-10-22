#include "coloBlendedViscosityRestrictionScheme.H"

#include "restrictionSchemes.H"
#include "addToRunTimeSelectionTable.H"

#ifdef NoRepository
#undef NoRepository
#endif

#include "colocated.H"
#include "twoPhaseModel.H"
#include "blendedViscosity.H"
#include "incompressibleTwoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

makeRestrictionSchemeNoTemplate(coloBlendedViscosity,lowerFaceScalar,colocated);

defineTemplateTypeNameAndDebugWithName
(
    blendedViscosityRestrictionScheme<colocated>,
    "blendedViscosity",
    0
);

coloBlendedViscosityRestrictionScheme::coloBlendedViscosityRestrictionScheme
(
    const fvMesh& fvMsh,
    Istream& is
)
:
    blendedViscosityRestrictionScheme<colocated>(fvMsh, is)
{
    const twoPhaseModel* modelPtr =
        &fvMsh.db().lookupObject<twoPhaseModel>("briscolaTwoPhaseDict");

    typedef blendedViscosity<incompressibleTwoPhaseModel<colocated>> modelType;

    if (dynamic_cast<const modelType*>(modelPtr))
    {
        this->C_ = dynamic_cast<const modelType*>(modelPtr)->C();
        this->mu1_ = dynamic_cast<const modelType*>(modelPtr)->mu1();
        this->mu2_ = dynamic_cast<const modelType*>(modelPtr)->mu2();
    }
    else
    {
        FatalErrorInFunction
            << "Two-phase model is of incompatible type"
            << endl << abort(FatalError);
    }
}

void coloBlendedViscosityRestrictionScheme::restrict
(
    meshDirection<lowerFaceScalar,colocated>& coarse,
    const meshDirection<lowerFaceScalar,colocated>& fine,
    const bool scale
)
{
    this->errorNoScaling(scale);

    const labelVector R(coarse.mshPart().R());

    coarse = Zero;

    const meshDirection<faceScalar,colocated>& faf =
        this->fvMsh().template
        metrics<colocated>().faceAreas()
        [fine.levelNum()][fine.directionNum()];

    forAllFaces(coarse, fd, i, j, k)
    {
        labelVector ijk(i,j,k);
        labelVector fijk(i*R.x(), j*R.y(), k*R.z());

        // Reconstruct the face alpha

        scalar alphaf = 0.0;
        scalar area = 0.0;

        labelVector o;
        for (o.x() = 0; o.x() < (fd == 0 ? 1 : R.x()); o.x()++)
        for (o.y() = 0; o.y() < (fd == 1 ? 1 : R.y()); o.y()++)
        for (o.z() = 0; o.z() < (fd == 2 ? 1 : R.z()); o.z()++)
        {
            alphaf += this->inv(fine(fijk+o)[fd])*faf(fijk+o)[fd*2];
            area += faf(fijk+o)[fd*2];
        }

        alphaf /= area;

        // Apply the blending function

        coarse(ijk)[fd] = this->blend(alphaf);
    }
}

}

}

}
