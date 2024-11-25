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
    avgColScalarFields_(),
    avgColVectorFields_(),
    avgColTensorFields_(),
    avgColSymmTensorFields_(),
    avgStagScalarFields_(),
    avgStagVectorFields_(),
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

    // Number of colocated/staggered fields needed
    label sizeColScalar = 0;
    label sizeColVector = 0;
    label sizeColTensor = 0;
    label sizeColSymmTensor = 0;
    label sizeStagScalar = 0;
    label sizeStagVector = 0;

    // Indices of field lists
    label indexColScalar = 0;
    label indexColVector = 0;
    label indexColTensor = 0;
    label indexColSymmTensor = 0;
    label indexStagScalar = 0;
    label indexStagVector = 0;

    forAll(fields_, i)
    {
        if (db.foundObject<colocatedScalarField>(fields_[i]))
        {
            sizeColScalar++;
        }
        else if (db.foundObject<colocatedVectorField>(fields_[i]))
        {
            sizeColVector++;
        }
        else if (db.foundObject<colocatedTensorField>(fields_[i]))
        {
            sizeColTensor++;
        }
        else if (db.foundObject<colocatedSymmTensorField>(fields_[i]))
        {
            sizeColSymmTensor++;
        }
        else if (db.foundObject<staggeredScalarField>(fields_[i]))
        {
            sizeStagScalar++;
        }
        else if (db.foundObject<staggeredVectorField>(fields_[i]))
        {
            sizeStagVector++;
        }
        else
        {
            WarningInFunction
                << "Field " << fields_[i] << " requested for sampling by "
                << this->name() << " but not found in registry." << endl;
        }
    }

    avgColScalarFields_.setSize(sizeColScalar);
    avgColVectorFields_.setSize(sizeColVector);
    avgColTensorFields_.setSize(sizeColTensor);
    avgColSymmTensorFields_.setSize(sizeColSymmTensor);

    avgStagScalarFields_.setSize(sizeStagScalar);
    avgStagVectorFields_.setSize(sizeStagVector);


    forAll(fields_, i)
    {
        if (db.foundObject<colocatedScalarField>(fields_[i]))
        {
            avgColScalarFields_.set
            (
                indexColScalar,
                new colocatedScalarField
                (
                    fields_[i]+"_avg",
                    fvMsh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    true
                )
            );

            indexColScalar++;
        }
        else if (db.foundObject<colocatedVectorField>(fields_[i]))
        {
            avgColVectorFields_.set
            (
                indexColVector,
                new colocatedVectorField
                (
                    fields_[i]+"_avg",
                    fvMsh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    true
                )
            );

            indexColVector++;
        }
        else if (db.foundObject<colocatedTensorField>(fields_[i]))
        {
            avgColTensorFields_.set
            (
                indexColTensor,
                new colocatedTensorField
                (
                    fields_[i]+"_avg",
                    fvMsh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    true
                )
            );

            indexColTensor++;
        }
        else if (db.foundObject<colocatedSymmTensorField>(fields_[i]))
        {
            avgColSymmTensorFields_.set
            (
                indexColSymmTensor,
                new colocatedSymmTensorField
                (
                    fields_[i]+"_avg",
                    fvMsh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    true
                )
            );

            indexColSymmTensor++;
        }
        else if (db.foundObject<staggeredScalarField>(fields_[i]))
        {
            avgStagScalarFields_.set
            (
                indexStagScalar,
                new staggeredScalarField
                (
                    fields_[i]+"_avg",
                    fvMsh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    true
                )
            );

            indexStagScalar++;
        }
        else if (db.foundObject<staggeredVectorField>(fields_[i]))
        {
            avgStagVectorFields_.set
            (
                indexStagVector,
                new staggeredVectorField
                (
                    fields_[i]+"_avg",
                    fvMsh_,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE,
                    true
                )
            );

            indexStagVector++;
        }
    }

    forAll(avgColScalarFields_, i)
    {
        avgColScalarFields_[i] = Zero;
    }

    forAll(avgColVectorFields_, i)
    {
        avgColVectorFields_[i] = Zero;
    }

    forAll(avgColTensorFields_, i)
    {
        avgColTensorFields_[i] = Zero;
    }

    forAll(avgColSymmTensorFields_, i)
    {
        avgColSymmTensorFields_[i] = Zero;
    }

    forAll(avgStagScalarFields_, i)
    {
        avgStagScalarFields_[i] = Zero;
    }

    forAll(avgStagVectorFields_, i)
    {
        avgStagVectorFields_[i] = Zero;
    }
}

timeAverage::~timeAverage()
{}

bool timeAverage::execute()
{
    if (runTime_.value() > startTime_)
    {
        const objectRegistry& db = fvMsh_.db();

        // Indices of field lists
        label indexColScalar = 0;
        label indexColVector = 0;
        label indexColTensor = 0;
        label indexColSymmTensor = 0;
        label indexStagScalar = 0;
        label indexStagVector = 0;

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
                    avgColScalarFields_[indexColScalar] = Zero;
                }

                forAllCells(csf[0],d,i,j,k)
                {
                    avgColScalarFields_[indexColScalar](0,d,i,j,k) *= avgT0;
                    avgColScalarFields_[indexColScalar](0,d,i,j,k) += csf(0,d,i,j,k)
                        * runTime_.deltaTValue();
                    avgColScalarFields_[indexColScalar](0,d,i,j,k) /= avgT;
                }

                indexColScalar++;
            }
            else if (db.foundObject<colocatedVectorField>(fields_[f]))
            {
                const colocatedVectorField& cvf
                    = db.lookupObject<colocatedVectorField>(fields_[f]);

                if (reset_)
                {
                    avgColVectorFields_[indexColVector] = Zero;
                }

                forAllCells(cvf[0],d,i,j,k)
                {
                    avgColVectorFields_[indexColVector](0,d,i,j,k) *= avgT0;
                    avgColVectorFields_[indexColVector](0,d,i,j,k) += cvf(0,d,i,j,k)
                        * runTime_.deltaTValue();
                    avgColVectorFields_[indexColVector](0,d,i,j,k) /= avgT;
                }

                indexColVector++;
            }
            else if (db.foundObject<colocatedTensorField>(fields_[f]))
            {
                const colocatedTensorField& ctf
                    = db.lookupObject<colocatedTensorField>(fields_[f]);

                if (reset_)
                {
                    avgColTensorFields_[indexColTensor] = Zero;
                }

                forAllCells(ctf[0],d,i,j,k)
                {
                    avgColTensorFields_[indexColTensor](0,d,i,j,k) *= avgT0;
                    avgColTensorFields_[indexColTensor](0,d,i,j,k) += ctf(0,d,i,j,k)
                        * runTime_.deltaTValue();
                    avgColTensorFields_[indexColTensor](0,d,i,j,k) /= avgT;
                }

                indexColTensor++;
            }
            else if (db.foundObject<colocatedSymmTensorField>(fields_[f]))
            {
                const colocatedSymmTensorField& cstf
                    = db.lookupObject<colocatedSymmTensorField>(fields_[f]);

                if (reset_)
                {
                    avgColSymmTensorFields_[indexColSymmTensor] = Zero;
                }

                forAllCells(cstf[0],d,i,j,k)
                {
                    avgColSymmTensorFields_[indexColSymmTensor](0,d,i,j,k) *= avgT0;
                    avgColSymmTensorFields_[indexColSymmTensor](0,d,i,j,k) += cstf(0,d,i,j,k)
                        * runTime_.deltaTValue();
                    avgColSymmTensorFields_[indexColSymmTensor](0,d,i,j,k) /= avgT;
                }

                indexColSymmTensor++;
            }
            else if (db.foundObject<staggeredScalarField>(fields_[f]))
            {
                const staggeredScalarField& ssf
                    = db.lookupObject<staggeredScalarField>(fields_[f]);

                if (reset_)
                {
                    avgStagScalarFields_[indexStagScalar] = Zero;
                }

                forAllCells(ssf[0],d,i,j,k)
                {
                    avgStagScalarFields_[indexStagScalar](0,d,i,j,k) *= avgT0;
                    avgStagScalarFields_[indexStagScalar](0,d,i,j,k)
                        += ssf(0,d,i,j,k) * runTime_.deltaTValue();
                    avgStagScalarFields_[indexStagScalar](0,d,i,j,k) /= avgT;
                }

                indexStagScalar++;
            }
            else if (db.foundObject<staggeredVectorField>(fields_[f]))
            {
                const staggeredVectorField& svf
                    = db.lookupObject<staggeredVectorField>(fields_[f]);

                if (reset_)
                {
                    avgStagVectorFields_[indexStagVector] = Zero;
                }

                forAllCells(svf[0],d,i,j,k)
                {
                    avgStagVectorFields_[indexStagVector](0,d,i,j,k) *= avgT0;
                    avgStagVectorFields_[indexStagVector](0,d,i,j,k)
                        += svf(0,d,i,j,k) * runTime_.deltaTValue();
                    avgStagVectorFields_[indexStagVector](0,d,i,j,k) /= avgT;
                }

                indexStagVector++;
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
