#include "budgetPrimitivesStaggered.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "exSchemesGradient.H"
#include "exSchemesReconstruction.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(budgetPrimitivesStaggered, 0);

addToRunTimeSelectionTable
(
    functionObject,
    budgetPrimitivesStaggered,
    dictionary
);

budgetPrimitivesStaggered::budgetPrimitivesStaggered
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict),
    fvMsh_(runTime.lookupObject<fvMesh>("briscolaMeshDict")),
    U2_
    (
        "U2",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    U3_
    (
        "U3",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    p2_
    (
        "p2",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    Rij_
    (
        "Rij",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    Rii_
    (
        "Rii",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    uiuiuk_
    (
        "uiuiuk",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    uiuiui_
    (
        "uiuiui",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    pui_
    (
        "pui",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    pdjui_
    (
        "pdjui",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    stagdjui_
    (
        "stagdjui",
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
        IOobject::AUTO_WRITE,
        true
    ),
    stagdjui2_
    (
        "stagdjui2",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        true
    ),
    djui2_
    (
        "djui2",
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
        IOobject::AUTO_WRITE,
        true
    ),
    eij_
    (
        "eij",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    )
{
    U2_        = Zero;
    U3_        = Zero;
    p2_        = Zero;
    Rij_       = Zero;
    Rii_       = Zero;
    uiuiuk_    = Zero;
    uiuiui_    = Zero;
    pui_       = Zero;
    pdjui_     = Zero;
    djui_      = Zero;
    stagdjui_  = Zero;
    djui2_     = Zero;
    stagdjui2_ = Zero;
    Sij_       = Zero;
    eij_       = Zero;
}

budgetPrimitivesStaggered::~budgetPrimitivesStaggered()
{}

bool budgetPrimitivesStaggered::execute()
{
    const objectRegistry& db = fvMsh_.db();

    const staggeredScalarField& U
        = db.lookupObject<staggeredScalarField>("U");

    const colocatedVectorField& Uc
        = db.lookupObject<colocatedVectorField>("Uc");

    const colocatedScalarField& p
        = db.lookupObject<colocatedScalarField>("p");

    forAllCells(U2_,l,d,i,j,k)
    {
        U2_(l,d,i,j,k) = Foam::sqr(U(l,d,i,j,k));
    }

    U2_.correct<bcsOfType<parallelBoundary>>();
    Rii_ = ex::reconstruct(U2_);

    forAllCells(U3_,l,d,i,j,k)
    {
        U3_(l,d,i,j,k) = Foam::pow3(U(l,d,i,j,k));
    }

    U3_.correct<bcsOfType<parallelBoundary>>();
    uiuiui_ = ex::reconstruct(U3_);

    forAllCells(p2_,l,d,i,j,k)
    {
        p2_(l,d,i,j,k) = Foam::sqr(p(l,d,i,j,k));
    }

    forAllCells(Rij_,l,d,i,j,k)
    {
        Rij_(l,d,i,j,k).xx() = Rii_(l,d,i,j,k).x();
        Rij_(l,d,i,j,k).xy() = Uc(l,d,i,j,k)[0]*Uc(l,d,i,j,k)[1];
        Rij_(l,d,i,j,k).xz() = Uc(l,d,i,j,k)[0]*Uc(l,d,i,j,k)[2];

        Rij_(l,d,i,j,k).yy() = Rii_(l,d,i,j,k).y();
        Rij_(l,d,i,j,k).yz() = Uc(l,d,i,j,k)[1]*Uc(l,d,i,j,k)[2];

        Rij_(l,d,i,j,k).zz() = Rii_(l,d,i,j,k).z();
    }

    forAllCells(uiuiuk_,l,d,i,j,k)
    {
        uiuiuk_(l,d,i,j,k).xx() = uiuiui_(l,d,i,j,k).x();
        uiuiuk_(l,d,i,j,k).xy() = Rij_(l,d,i,j,k).xx()*Uc(l,d,i,j,k)[1];
        uiuiuk_(l,d,i,j,k).xz() = Rij_(l,d,i,j,k).xx()*Uc(l,d,i,j,k)[2];

        uiuiuk_(l,d,i,j,k).yx() = Rij_(l,d,i,j,k).yy()*Uc(l,d,i,j,k)[0];
        uiuiuk_(l,d,i,j,k).yy() = uiuiui_(l,d,i,j,k).y();
        uiuiuk_(l,d,i,j,k).yz() = Rij_(l,d,i,j,k).yy()*Uc(l,d,i,j,k)[2];

        uiuiuk_(l,d,i,j,k).zx() = Rij_(l,d,i,j,k).zz()*Uc(l,d,i,j,k)[0];
        uiuiuk_(l,d,i,j,k).zy() = Rij_(l,d,i,j,k).zz()*Uc(l,d,i,j,k)[1];
        uiuiuk_(l,d,i,j,k).zz() = uiuiui_(l,d,i,j,k).z();
    }

    forAllCells(pui_,l,d,i,j,k)
    {
        pui_(l,d,i,j,k)[0] = p(l,d,i,j,k)*Uc(l,d,i,j,k)[0];
        pui_(l,d,i,j,k)[1] = p(l,d,i,j,k)*Uc(l,d,i,j,k)[1];
        pui_(l,d,i,j,k)[2] = p(l,d,i,j,k)*Uc(l,d,i,j,k)[2];
    }

    stagdjui_ = ex::grad(U);
    stagdjui_.correct<bcsOfType<parallelBoundary>>();
    djui_ = ex::reconstruct(stagdjui_);

    forAllCells(pdjui_,l,d,i,j,k)
    {
        pdjui_(l,d,i,j,k) = p(l,d,i,j,k)*djui_(l,d,i,j,k);
    }

    forAllCells(stagdjui2_, l, d, i, j, k)
    {
        stagdjui2_(l,d,i,j,k) = Foam::sqr(stagdjui_(l,d,i,j,k)[d]);
    }

    stagdjui2_.correct<bcsOfType<parallelBoundary>>();
    djui2_ = ex::reconstruct(stagdjui2_);

    forAllCells(Sij_,l,d,i,j,k)
    {
        Sij_(l,d,i,j,k) = 0.5*(djui_(l,d,i,j,k)+djui_(l,d,i,j,k).T());
    }

    forAllCells(eij_,l,d,i,j,k)
    {
        eij_(l,d,i,j,k).xx() = djui2_(l,d,i,j,k).x();
        eij_(l,d,i,j,k).xy() = Sij_(l,d,i,j,k).xy()*djui_(l,d,i,j,k).xy();
        eij_(l,d,i,j,k).xz() = Sij_(l,d,i,j,k).xz()*djui_(l,d,i,j,k).xz();

        eij_(l,d,i,j,k).yx() = Sij_(l,d,i,j,k).yx()*djui_(l,d,i,j,k).yx();
        eij_(l,d,i,j,k).yy() = djui2_(l,d,i,j,k).y();
        eij_(l,d,i,j,k).yz() = Sij_(l,d,i,j,k).yz()*djui_(l,d,i,j,k).yz();

        eij_(l,d,i,j,k).zx() = Sij_(l,d,i,j,k).zx()*djui_(l,d,i,j,k).zx();
        eij_(l,d,i,j,k).zy() = Sij_(l,d,i,j,k).zy()*djui_(l,d,i,j,k).zy();
        eij_(l,d,i,j,k).zz() = djui2_(l,d,i,j,k).z();
    }

    return true;
}

bool budgetPrimitivesStaggered::write()
{
    return true;
}

bool budgetPrimitivesStaggered::end()
{
    return true;
}

}

}

}

}
