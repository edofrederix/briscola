#include "phaseInversion.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "exSchemes.H"
#include <fstream>
#include "PstreamReduceOps.H"

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
    briscolaFunctionObject(name, runTime, dict),
    fvMsh_(runTime.lookupObject<fvMesh>("briscolaMeshDict")),
    twoPhaseDict_
    (
        IOobject
        (
            "briscolaTwoPhaseDict",
            runTime.system(),
            runTime,
            Foam::IOobject::MUST_READ,
            Foam::IOobject::NO_WRITE
        )
    ),
    rho1_(twoPhaseDict_.lookup<scalar>("rho1")),
    rho2_(twoPhaseDict_.lookup<scalar>("rho2")),
    g_(twoPhaseDict_.lookup<vector>("g")),
    magSqrU_
    (
        "magSqrU",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    y_
    (
        "y",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
{
    read(dict);

    if (Pstream::master())
    {
        filePtr_.reset(new OFstream("data.txt"));
    }

    const colocatedVectorField& cc =
        fvMsh_.metrics<colocated>().cellCenters();

    y_ = Zero;

    forAllCells(y_,i,j,k)
    {
        y_(i,j,k) = cc(i,j,k).y();
    }
}

phaseInversion::~phaseInversion()
{}

bool phaseInversion::read(const dictionary& dict)
{
    if (runTime_.time().value() == 0.0)
    {
        colocatedScalarField& alpha =
            runTime_.lookupObjectRef<colocatedScalarField>("alpha");

        alpha = Zero;

        const colocatedVectorField& cc =
            fvMsh_.metrics<colocated>().cellCenters();

        const double h = 0.0; // cube offset

        forAllCells(alpha, i, j, k)
            if (
                    cc(i,j,k).x() > h && cc(i,j,k).x() < (0.5 + h) &&
                    cc(i,j,k).y() > h && cc(i,j,k).y() < (0.5 + h) &&
                    cc(i,j,k).z() > h && cc(i,j,k).z() < (0.5 + h)
                )
                {
                    alpha(i,j,k) = 1.0;
                }


        alpha.correctBoundaryConditions();
    }

    return true;
}

bool phaseInversion::execute()
{
    const colocatedScalarDirection& alpha =
        runTime_.lookupObjectRef<colocatedScalarField>("alpha")[0][0];

    const colocatedVectorDirection& Uc =
        runTime_.lookupObjectRef<colocatedVectorField>("Uc")[0][0];

    const colocatedScalarDirection& cv =
        fvMsh_.metrics<colocated>().cellVolumes()[0][0];

    forAllCells(Uc, i, j, k)
    {
        magSqrU_(i,j,k) = magSqr(Uc(i,j,k));
    }

    const scalar ke1(gSum(cv*alpha*0.5*rho1_*magSqrU_[0][0]));
    const scalar ke2(gSum(cv*(1.0-alpha)*0.5*rho2_*magSqrU_[0][0]));

    const scalar pe1(gSum(cv*alpha*rho1_*mag(g_)*y_[0][0]));
    const scalar pe2(gSum(cv*(1.0-alpha)*rho2_*mag(g_)*y_[0][0]));

    // Normalizations

    const scalar tstar = runTime_.time().value() * 0.245;

    const scalar ke1star = ke1 / (rho1_ * Foam::sqr(0.245) / 16.0);
    const scalar ke2star = ke2 / (rho2_ * Foam::sqr(0.245) / 16.0);

    const scalar pe1star = pe1 / (15.0/128.0 * rho1_ * Foam::mag(g_));
    const scalar pe2star = pe2 / (49.0/128.0 * rho2_ * Foam::mag(g_));

    if (Pstream::master())
        filePtr_() << tstar
                << " " << ke1star << " " << ke2star
                << " " << pe1star << " " << pe2star
                << nl;

    return true;
}

}

}

}

}
