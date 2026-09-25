#include "budgetPrimitivesStaggeredTwoPhase.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

#undef NoRepository
#include "incompressibleTwoPhaseModel.H"
#define NoRepository

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(budgetPrimitivesStaggeredTwoPhase, 0);

addToRunTimeSelectionTable
(
    functionObject,
    budgetPrimitivesStaggeredTwoPhase,
    dictionary
);

budgetPrimitivesStaggeredTwoPhase::budgetPrimitivesStaggeredTwoPhase
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict),
    samplingSteps_
    (
        dict.lookupOrDefault<label>("samplingSteps", 1)
    ),
    fvMsh_(runTime.lookupObject<fvMesh>("briscolaMeshDict")),
    tpm_
    (
        fvMsh_.db().lookupObject<twoPhaseModel>("briscolaTwoPhaseDict")
    ),
    Uc_(fvMsh_.db().lookupObject<colocatedVectorField>("Uc")),
    p_(fvMsh_.db().lookupObject<colocatedScalarField>("p")),
    f1_
    (
        "f1",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    f2_(tpm_.alpha()),

    rho1_(tpm_.cast<const incompressibleTwoPhaseModel<staggered>>().rho1()),
    rho2_(tpm_.cast<const incompressibleTwoPhaseModel<staggered>>().rho2()),

    mu1_(tpm_.dict().lookup<scalar>("mu1")),
    mu2_(tpm_.dict().lookup<scalar>("mu2")),
    rho_
    (
        "rho",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    mu_
    (
        "mu",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    si_
    (
        "si",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    djui_
    (
        "djui",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    Sij_
    (
        "Sij",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    bp1_
    (
        *this,
        f1_,
        1
    ),
    bp2_
    (
        *this,
        f2_,
        2
    )
{
    f1_             = Zero;
    si_             = Zero;
    djui_           = Zero;
    Sij_            = Zero;
}

budgetPrimitivesStaggeredTwoPhase::~budgetPrimitivesStaggeredTwoPhase()
{}

bool budgetPrimitivesStaggeredTwoPhase::execute()
{
    if (runTime_.timeIndex() % samplingSteps_ == 0)
    {
        f1_ = 1.0 - f2_;

        djui_ = ex::grad(Uc_);
        djui_.correctBoundaryConditions();

        Sij_ = 0.5*(djui_ + T(djui_));

        rho_ = f1_ * rho1_ + f2_ * rho2_;

        mu_ = f1_ * mu1_ + f2_ * mu2_;

        si_ = ex::reconstruct(tpm_.alphavf().surfaceTension());

        bp1_.correct();
        bp2_.correct();
    }
    return true;
}

bool budgetPrimitivesStaggeredTwoPhase::write()
{
    return true;
}

bool budgetPrimitivesStaggeredTwoPhase::end()
{
    return true;
}

}

}

}

}
