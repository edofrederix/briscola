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
    name_(name),
    timeAveragedFields_(fields_.size()),
    startTime_
    (
        Foam::max
        (
            dict.lookupOrDefault<scalar>("startTime", 0.0),
            runTime.startTime().value()
        )
    )
{
    init();
}

void timeAverage::init()
{
    const objectRegistry& db = fvMsh_.db();

    scalar size = 0;
    scalar index = 0;

    forAll(fields_, i)
    {
        if (db.foundObject<colocatedScalarField>(fields_[i]))
        {
            size += 1;
        }
        else if (db.foundObject<colocatedVectorField>(fields_[i]))
        {
            size += 3;
        }
        else
        {
            WarningInFunction
                << "Field " << fields_[i] << " requested for sampling by "
                << name_ << " but not found in registry." << endl;
        }
    }

    timeAveragedFields_.setSize(size);

    forAll(fields_, i)
    {
        if (db.foundObject<colocatedScalarField>(fields_[i]))
        {
            timeAveragedFields_.set
            (
                index,
                new colocatedScalarField
                (
                    fields_[i]+"_avg",
                    fvMsh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    true
                )
            );

            index++;
        }
        else if (db.foundObject<colocatedVectorField>(fields_[i]))
        {
            for (int d = 0; d < 3; d++)
            {
                word dir(Foam::name(d));

                timeAveragedFields_.set
                (
                    index,
                    new colocatedScalarField
                    (
                        fields_[i]+dir+"_avg",
                        fvMsh_,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE,
                        true
                    )
                );

                index++;
            }
        }
    }

    forAll(timeAveragedFields_, i)
    {
        timeAveragedFields_[i] = Zero;
    }
}

timeAverage::~timeAverage()
{}

bool timeAverage::execute()
{
    if (runTime_.value() > startTime_)
    {
        const objectRegistry& db = fvMsh_.db();

        scalar index = 0;

        scalar avgT  = runTime_.value() - startTime_;
        scalar avgT0 = avgT - runTime_.deltaTValue();


        forAll(fields_, f)
        {
            if (db.foundObject<colocatedScalarField>(fields_[f]))
            {
                const colocatedScalarField& csf
                    = db.lookupObject<colocatedScalarField>(fields_[f]);

                forAllLevels(csf,l,d,i,j,k)
                {
                    timeAveragedFields_[index](l,d,i,j,k) *= avgT0;
                    timeAveragedFields_[index](l,d,i,j,k) += csf(l,d,i,j,k)
                        * runTime_.deltaTValue();
                    timeAveragedFields_[index](l,d,i,j,k) /= avgT;
                }

                index++;
            }
            else if (db.foundObject<colocatedVectorField>(fields_[f]))
            {
                const colocatedVectorField& cvf
                    = db.lookupObject<colocatedVectorField>(fields_[f]);

                for (int dir = 0; dir < 3; dir++)
                {
                    forAllLevels(cvf,l,d,i,j,k)
                    {
                        timeAveragedFields_[index](l,d,i,j,k) *= avgT0;
                        timeAveragedFields_[index](l,d,i,j,k)
                            += cvf(l,d,i,j,k)[dir] * runTime_.deltaTValue();
                        timeAveragedFields_[index](l,d,i,j,k) /= avgT;
                    }

                    index++;
                }
            }
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
