#include "immersedBoundaryForceStaggered.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "IOdictionary.H"
#include "exSchemes.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(immersedBoundaryForceStaggered, 0);

addToRunTimeSelectionTable
(
    functionObject,
    immersedBoundaryForceStaggered,
    dictionary
);

immersedBoundaryForceStaggered::immersedBoundaryForceStaggered
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict),
    fvMsh_(runTime.lookupObject<fvMesh>("briscolaMeshDict")),
    volForce_
    (
        "volForce",
        fvMsh_,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE,
        true
    ),
    name_(word(dict.lookup("name"))),
    method_(word(dict.lookup("method"))),
    exSourceForce_(Zero)
{
    bool found = false;

    forAll(fvMsh_.IB<staggered>(), ib)
    {
        if (fvMsh_.IB<staggered>()[ib].name() == name_)
        {
            number_ = ib;
            found = true;
        }
    }

    if (!found)
    {
        FatalError
            << "Immersed boundary " << name_
            << " requested but not found."
            << endl;
        FatalError.exit();
    }

    if
    (
        !
        (
            method_ == "penalization"
            || method_ == "Mittal"
            || method_ == "Vreman"
            || method_ == "Fadlun"
        )
    )
    {
        FatalError
            << "Unknown immersed boundary method: "
            << method_
            << endl;
        FatalError.exit();
    }

    computeExplicitSourceForce();
}

immersedBoundaryForceStaggered::~immersedBoundaryForceStaggered()
{}

bool immersedBoundaryForceStaggered::execute()
{
    const IOdictionary& solverDict
        = runTime_.lookupObject<IOdictionary>("briscolaStaggeredDict");

    const scalar nu = readScalar(solverDict.lookup("nu"));

    const objectRegistry& db = fvMsh_.db();

    const colocatedScalarField& p
        = db.lookupObject<colocatedScalarField>("p");

    const staggeredScalarField& U
        = db.lookupObject<staggeredScalarField>("U");

    staggeredLowerFaceScalarField phi = ex::faceFlux(U);

    const staggeredScalarField& cv
        = fvMsh_.metrics<staggered>().cellVolumes();

    volForce_ = Zero;
    volForce_ +=
        - ex::laplacian(nu,U)
        + ex::div(phi,U)
        + ex::stagGrad(p);

    if (method_ == "Mittal" || method_ == "Vreman")
    {
        volForce_ *= fvMsh_.IB<staggered>()[number_].ghostMask();
    }
    else if (method_ == "penalization")
    {
        volForce_ *= fvMsh_.IB<staggered>()[number_].mask();
    }
    else if (method_ == "Fadlun")
    {
        volForce_ *= fvMsh_.IB<staggered>()[number_].wallAdjMask();
    }
    else
    {
        FatalError
            << "Unknown immersed boundary method: "
            << method_
            << endl;
        FatalError.exit();
    }

    vector totalForce = Zero;

    forAllCells(volForce_[0],d,i,j,k)
    {
        totalForce[d] += cv(0,d,i,j,k)*volForce_(0,d,i,j,k);
    }

    totalForce -= exSourceForce_;

    // Gather all total forces on master proc

    autoPtr<List<vector>> totalForcesList;

    if (Pstream::master())
    {
        totalForcesList.reset(new List<vector>(Pstream::nProcs(), Zero));
    }

    // Receive sizes and displacements
    labelList rs(Pstream::nProcs(), sizeof(vector));
    labelList rd(Pstream::nProcs(), Zero);

    forAll(rd, p)
    {
        rd[p] = p*rs[p];
    }

    UPstream::gather
    (
        reinterpret_cast<char*>(&totalForce),
        sizeof(vector),
        reinterpret_cast<char*>
        (
            Pstream::master()
        ? totalForcesList->begin()
        : nullptr
        ),
        rs,
        rd,
        UPstream::worldComm
    );

    if (Pstream::master())
    {
        vector globalTotalForce = Zero;

        forAll(totalForcesList(), f)
        {
            globalTotalForce += totalForcesList()[f];
        }
        Info << "IBM force on " << name_ << " = " << globalTotalForce << endl;
    }

    return true;
}

bool immersedBoundaryForceStaggered::write()
{
    return true;
}

bool immersedBoundaryForceStaggered::end()
{
    return true;
}

void immersedBoundaryForceStaggered::computeExplicitSourceForce()
{
    const IOdictionary& solverDict
        = runTime_.lookupObject<IOdictionary>("briscolaStaggeredDict");

    const vector exSource = vector(solverDict.lookup("source"));

    const staggeredScalarField& cv
        = fvMsh_.metrics<staggered>().cellVolumes();

    if (mag(exSource) > 1e-10)
    {
        forAllCells(fvMsh_.IB<staggered>()[number_].mask()[0],d,i,j,k)
        {
            exSourceForce_[d] += exSource[d]
                * cv(0,d,i,j,k)
                * fvMsh_.IB<staggered>()[number_].mask()(0,d,i,j,k);

        }
        if (method_ == "Fadlun")
        {
            forAllCells(fvMsh_.IB<staggered>()[number_].wallAdjMask()[0],d,i,j,k)
            {
                exSourceForce_[d] += exSource[d]
                    * cv(0,d,i,j,k)
                    * fvMsh_.IB<staggered>()[number_].wallAdjMask()(0,d,i,j,k);

            }
        }
    }
}

}

}

}

}
