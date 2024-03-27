#include "budgetPrimitives.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(budgetPrimitives, 0);

addToRunTimeSelectionTable
(
    functionObject,
    budgetPrimitives,
    dictionary
);

budgetPrimitives::budgetPrimitives
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict),
    fvMsh_(runTime.lookupObject<fvMesh>("briscolaMeshDict")),
    uName_
    (
        fvMsh_.db().foundObject<colocatedVectorField>("Uc") ?
        "Uc" : "U"
    ),
    u2_
    (
        "u2",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    v2_
    (
        "v2",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    w2_
    (
        "w2",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    p2_
    (
        "p2",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    )
{
    u2_ = Zero;
    v2_ = Zero;
    w2_ = Zero;
    p2_ = Zero;

    if (!fvMsh_.db().foundObject<colocatedVectorField>(uName_))
    {
        FatalError
            << "Velocity field " << uName_ << " not found."
            << endl << abort(FatalError);
    }
}

budgetPrimitives::~budgetPrimitives()
{}

bool budgetPrimitives::execute()
{
    const objectRegistry& db = fvMsh_.db();

    const colocatedVectorField& U
        = db.lookupObject<colocatedVectorField>(uName_);

    const colocatedScalarField& p
        = db.lookupObject<colocatedScalarField>("p");

    forAllCells(u2_,l,d,i,j,k)
    {
        u2_(l,d,i,j,k) = Foam::sqr(U(l,d,i,j,k)[0]);
        v2_(l,d,i,j,k) = Foam::sqr(U(l,d,i,j,k)[1]);
        w2_(l,d,i,j,k) = Foam::sqr(U(l,d,i,j,k)[2]);
        p2_(l,d,i,j,k) = Foam::sqr(p(l,d,i,j,k));
    }

    return true;
}

bool budgetPrimitives::write()
{
    return true;
}

bool budgetPrimitives::end()
{
    return true;
}

}

}

}

}
