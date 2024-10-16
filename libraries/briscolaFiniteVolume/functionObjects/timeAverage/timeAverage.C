#include "timeAverage.H"
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

defineTypeNameAndDebug(timeAverage, 0);

addToRunTimeSelectionTable
(
    functionObject,
    timeAverage,
    dictionary
);

timeAverage::timeAverage
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    briscolaFunctionObject(name, runTime, dict),
    fvMsh_(runTime.lookupObject<fvMesh>("briscolaMeshDict")),
    fields_(dict.lookup("fields")),
    timeAveragedColFields_(),
    timeAveragedStagFields_(),
    startTime_
    (
        Foam::max
        (
            dict.lookupOrDefault<scalar>("startTime", 0.0),
            runTime.startTime().value()
        )
    ),
    intervalAveraging_(dict.lookupOrDefault("intervalAveraging", false)),
    intervalStart_(startTime_),
    reset_(false)
{
    init();
}

void timeAverage::init()
{
    const objectRegistry& db = fvMsh_.db();

    // Number of colocated/staggered scalar fields needed
    scalar sizeCol = 0;
    scalar sizeStag = 0;

    // Indices of field lists
    scalar indexCol = 0;
    scalar indexStag = 0;

    forAll(fields_, i)
    {
        if (db.foundObject<colocatedScalarField>(fields_[i]))
        {
            sizeCol += 1;
        }
        else if (db.foundObject<colocatedVectorField>(fields_[i]))
        {
            sizeCol += 3;
        }
        else if (db.foundObject<colocatedSymmTensorField>(fields_[i]))
        {
            sizeCol += 6;
        }
        else if (db.foundObject<colocatedTensorField>(fields_[i]))
        {
            sizeCol += 9;
        }
        else if (db.foundObject<staggeredScalarField>(fields_[i]))
        {
            sizeStag += 1;
        }
        else
        {
            WarningInFunction
                << "Field " << fields_[i] << " requested for sampling by "
                << this->name() << " but not found in registry." << endl;
        }
    }

    timeAveragedColFields_.setSize(sizeCol);
    timeAveragedStagFields_.setSize(sizeStag);

    forAll(fields_, i)
    {
        if (db.foundObject<colocatedScalarField>(fields_[i]))
        {
            timeAveragedColFields_.set
            (
                indexCol,
                new colocatedScalarField
                (
                    fields_[i]+"_avg",
                    fvMsh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    true
                )
            );

            indexCol++;
        }
        else if (db.foundObject<colocatedVectorField>(fields_[i]))
        {
            for (int c = 0; c < 3; c++)
            {
                wordList cmpts({"X", "Y", "Z"});

                timeAveragedColFields_.set
                (
                    indexCol,
                    new colocatedScalarField
                    (
                        fields_[i]+"_"+cmpts[c]+"_avg",
                        fvMsh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        true
                    )
                );

                indexCol++;
            }
        }
        else if (db.foundObject<colocatedSymmTensorField>(fields_[i]))
        {
            wordList cmpts({"XX", "XY", "XZ", "YY", "YZ", "ZZ"});

            for (int c = 0; c < 6; c++)
            {
                timeAveragedColFields_.set
                (
                    indexCol,
                    new colocatedScalarField
                    (
                        fields_[i]+"_"+cmpts[c]+"_avg",
                        fvMsh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        true
                    )
                );

                indexCol++;
            }
        }
        else if (db.foundObject<colocatedTensorField>(fields_[i]))
        {
            wordList cmpts
            (
                {"XX", "XY", "XZ", "YX", "YY", "YZ", "ZX", "ZY", "ZZ"}
            );

            for (int c = 0; c < 9; c++)
            {
                timeAveragedColFields_.set
                (
                    indexCol,
                    new colocatedScalarField
                    (
                        fields_[i]+"_"+cmpts[c]+"_avg",
                        fvMsh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        true
                    )
                );

                indexCol++;
            }
        }
        else if (db.foundObject<staggeredScalarField>(fields_[i]))
        {
            timeAveragedStagFields_.set
            (
                indexStag,
                new staggeredScalarField
                (
                    fields_[i]+"_avg",
                    fvMsh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    true
                )
            );

            indexStag++;
        }
    }

    forAll(timeAveragedColFields_, i)
    {
        timeAveragedColFields_[i] = Zero;
    }

    forAll(timeAveragedStagFields_, i)
    {
        timeAveragedStagFields_[i] = Zero;
    }
}

timeAverage::~timeAverage()
{}

bool timeAverage::execute()
{
    if (runTime_.value() > startTime_)
    {
        const objectRegistry& db = fvMsh_.db();

        scalar indexCol = 0;
        scalar indexStag = 0;

        scalar avgT  = runTime_.value() - intervalStart_;
        scalar avgT0 = avgT - runTime_.deltaTValue();

        forAll(fields_, f)
        {
            if (db.foundObject<colocatedScalarField>(fields_[f]))
            {
                const colocatedScalarField& csf
                    = db.lookupObject<colocatedScalarField>(fields_[f]);

                if (reset_)
                {
                    timeAveragedColFields_[indexCol] = Zero;
                }

                forAllCells(csf,l,d,i,j,k)
                {
                    timeAveragedColFields_[indexCol](l,d,i,j,k) *= avgT0;
                    timeAveragedColFields_[indexCol](l,d,i,j,k) += csf(l,d,i,j,k)
                        * runTime_.deltaTValue();
                    timeAveragedColFields_[indexCol](l,d,i,j,k) /= avgT;
                }

                indexCol++;
            }
            else if (db.foundObject<colocatedVectorField>(fields_[f]))
            {
                const colocatedVectorField& cvf
                    = db.lookupObject<colocatedVectorField>(fields_[f]);

                for (int c = 0; c < 3; c++)
                {
                    if (reset_)
                    {
                        timeAveragedColFields_[indexCol] = Zero;
                    }

                    forAllCells(cvf,l,d,i,j,k)
                    {
                        timeAveragedColFields_[indexCol](l,d,i,j,k) *= avgT0;
                        timeAveragedColFields_[indexCol](l,d,i,j,k)
                            += cvf(l,d,i,j,k)[c] * runTime_.deltaTValue();
                        timeAveragedColFields_[indexCol](l,d,i,j,k) /= avgT;
                    }

                    indexCol++;
                }
            }
            else if (db.foundObject<colocatedSymmTensorField>(fields_[f]))
            {
                const colocatedSymmTensorField& cstf
                    = db.lookupObject<colocatedSymmTensorField>(fields_[f]);

                for (int c = 0; c < 6; c++)
                {
                    if (reset_)
                    {
                        timeAveragedColFields_[indexCol] = Zero;
                    }

                    forAllCells(cstf,l,d,i,j,k)
                    {
                        timeAveragedColFields_[indexCol](l,d,i,j,k) *= avgT0;
                        timeAveragedColFields_[indexCol](l,d,i,j,k)
                            += cstf(l,d,i,j,k)[c] * runTime_.deltaTValue();
                        timeAveragedColFields_[indexCol](l,d,i,j,k) /= avgT;
                    }

                    indexCol++;
                }
            }
            else if (db.foundObject<colocatedTensorField>(fields_[f]))
            {
                const colocatedTensorField& ctf
                    = db.lookupObject<colocatedTensorField>(fields_[f]);

                for (int c = 0; c < 9; c++)
                {
                    if (reset_)
                    {
                        timeAveragedColFields_[indexCol] = Zero;
                    }

                    forAllCells(ctf,l,d,i,j,k)
                    {
                        timeAveragedColFields_[indexCol](l,d,i,j,k) *= avgT0;
                        timeAveragedColFields_[indexCol](l,d,i,j,k)
                            += ctf(l,d,i,j,k)[c] * runTime_.deltaTValue();
                        timeAveragedColFields_[indexCol](l,d,i,j,k) /= avgT;
                    }

                    indexCol++;
                }
            }
            else if (db.foundObject<staggeredScalarField>(fields_[f]))
            {
                const staggeredScalarField& ssf
                    = db.lookupObject<staggeredScalarField>(fields_[f]);

                if (reset_)
                {
                    timeAveragedStagFields_[indexStag] = Zero;
                }

                forAllCells(ssf,l,d,i,j,k)
                {
                    timeAveragedStagFields_[indexStag](l,d,i,j,k) *= avgT0;
                    timeAveragedStagFields_[indexStag](l,d,i,j,k)
                        += ssf(l,d,i,j,k) * runTime_.deltaTValue();
                    timeAveragedStagFields_[indexStag](l,d,i,j,k) /= avgT;
                }

                indexStag++;
            }
        }

        reset_ = false;

        if (runTime_.writeTime() && intervalAveraging_)
        {
            intervalStart_ = runTime_.value();
            reset_ = true;
        }
    }

    return true;
}

bool timeAverage::write()
{
    return true;
}

bool timeAverage::end()
{
    return true;
}

}

}

}

}
