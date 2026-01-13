#include "phaseInversion.H"
#include "addToRunTimeSelectionTable.H"
#include "incompressibleTwoPhaseModel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(phaseInversion, 0);

addToRunTimeSelectionTable
(
    functionObject,
    phaseInversion,
    dictionary
);

phaseInversion::phaseInversion
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict)
{
    read(dict);

    if (Pstream::master())
        filePtr_.reset(new OFstream("data.txt"));

    // Read two-phase properties

    const twoPhaseModel& model =
        runTime.lookupObject<twoPhaseModel>("briscolaTwoPhaseDict");

    g_ = Foam::mag(model.g());

    if (model.castable<incompressibleTwoPhaseModel<colocated>>())
    {
        const incompressibleTwoPhaseModel<colocated>& iModel =
            model.cast<incompressibleTwoPhaseModel<colocated>>();

        rho1_ = iModel.rho1();
        rho2_ = iModel.rho2();
    }
    else
    {
        const incompressibleTwoPhaseModel<staggered>& iModel =
            model.cast<incompressibleTwoPhaseModel<staggered>>();

        rho1_ = iModel.rho1();
        rho2_ = iModel.rho2();
    }
}

bool phaseInversion::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        alpha = Zero;

        const fvMesh& fvMsh = alpha.fvMsh();

        const colocatedVectorField& cc =
            fvMsh.metrics<colocated>().cellCenters();

        forAllCells(alpha, i, j, k)
        {
            if
            (
                cc(i,j,k).x() > 0 && cc(i,j,k).x() < 0.5
             && cc(i,j,k).y() > 0 && cc(i,j,k).y() < 0.5
             && cc(i,j,k).z() > 0 && cc(i,j,k).z() < 0.5
            )
            {
                alpha(i,j,k) = 1.0;
            }
        }

        alpha.correctBoundaryConditions();
    }

    return true;
}

bool phaseInversion::execute()
{
    const colocatedScalarField& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha");

    const colocatedVectorField& U =
        runTime_.foundObject<colocatedVectorField>("U")
      ? runTime_.lookupObject<colocatedVectorField>("U")
      : runTime_.lookupObject<colocatedVectorField>("Uc");

    const fvMesh& fvMsh = U.fvMsh();

    const colocatedScalarField& cv =
        fvMsh.metrics<colocated>().cellVolumes();

    const colocatedVectorField& cc =
        fvMsh.metrics<colocated>().cellCenters();

    const colocatedScalarField magSqrU(mag(U & U));

    const colocatedScalarField y(cc.component(vector::Y));

    const scalar ke1(gSum(cv*alpha*0.5*rho1_*magSqrU)[0]);
    const scalar ke2(gSum(cv*(1.0-alpha)*0.5*rho2_*magSqrU)[0]);

    const scalar pe1(gSum(cv*alpha*rho1_*g_*y)[0]);
    const scalar pe2(gSum(cv*(1.0-alpha)*rho2_*g_*y)[0]);

    // Normalizations

    const scalar tStar = runTime_.time().value()*0.245;

    const scalar ke1star = ke1/(rho1_*Foam::sqr(0.245)/16.0);
    const scalar ke2star = ke2/(rho2_*Foam::sqr(0.245)/16.0);

    const scalar pe1star = pe1/(15.0/128.0*rho1_*g_);
    const scalar pe2star = pe2/(49.0/128.0*rho2_*g_);

    if (Pstream::master())
        filePtr_()
            << tStar
            << " " << ke1star << " " << ke2star
            << " " << pe1star << " " << pe2star
            << nl;

    return true;
}

}

}

}

}
