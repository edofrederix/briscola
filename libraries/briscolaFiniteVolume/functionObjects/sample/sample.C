#include "sample.H"
#include "Time.H"
#include "OSspecific.H"
#include "OFstream.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

void sample::init()
{
    interpPtr_.reset
    (
        pointInterpolator<colocated>::New
        (
            fvMsh_,
            this->points(),
            dict_.lookupOrDefault<word>("interpolator", "linear")
        ).ptr()
    );

    if (!interpPtr_->allPointsFound())
        WarningInFunction
            << "Not all points were found for " << name_ << " sample."
            << endl;



    const objectRegistry& db = fvMsh_.db();

    scalar size = 0;

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

    timeAveragedData_.setSize(timeAverage_*size);

    forAll(timeAveragedData_, i)
    {
        timeAveragedData_.set
        (
            i,
            new scalarList
            (
                interpPtr_->points().size(),
                Zero
            )
        );
    }
}

void sample::appendData
(
    const word fieldName,
    PtrList<scalarList>& data,
    wordList& headers
)
{
    const objectRegistry& db = fvMsh_.db();

    if (db.foundObject<colocatedScalarField>(fieldName))
    {
        this->appendData
        (
            db.lookupObject<colocatedScalarField>(fieldName),
            data,
            headers
        );
    }
    else if (db.foundObject<colocatedVectorField>(fieldName))
    {
        this->appendData
        (
            db.lookupObject<colocatedVectorField>(fieldName),
            data,
            headers
        );
    }
    else if (db.foundObject<colocatedTensorField>(fieldName))
    {
        this->appendData
        (
            db.lookupObject<colocatedTensorField>(fieldName),
            data,
            headers
        );
    }
    else if (db.foundObject<colocatedSphericalTensorField>(fieldName))
    {
        this->appendData
        (
            db.lookupObject<colocatedSphericalTensorField>(fieldName),
            data,
            headers
        );
    }
    else if (db.foundObject<colocatedSymmTensorField>(fieldName))
    {
        this->appendData
        (
            db.lookupObject<colocatedSymmTensorField>(fieldName),
            data,
            headers
        );
    }
    else if (db.foundObject<colocatedDiagTensorField>(fieldName))
    {
        this->appendData
        (
            db.lookupObject<colocatedDiagTensorField>(fieldName),
            data,
            headers
        );
    }
    else if (db.foundObject<colocatedFaceScalarField>(fieldName))
    {
        this->appendData
        (
            db.lookupObject<colocatedFaceScalarField>(fieldName),
            data,
            headers
        );
    }
    else if (db.foundObject<colocatedFaceVectorField>(fieldName))
    {
        this->appendData
        (
            db.lookupObject<colocatedFaceVectorField>(fieldName),
            data,
            headers
        );
    }
    else
    {
        WarningInFunction
            << "Field " << fieldName << " requested for sampling by "
            << name_ << " but not found in registry." << endl;
    }
}

sample::sample
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
    timeAverage_(dict.lookupOrDefault<Switch>("timeAverage", false)),
    timeAveragedData_(timeAverage_*fields_.size()),
    startTime_
    (
        Foam::max
        (
            dict.lookupOrDefault<scalar>("startTime", 0.0),
            runTime.startTime().value()
        )
    ),
    writeTime_(Foam::name(startTime_))
{
    init();
}

sample::~sample()
{}

bool sample::execute()
{
    return true;
}

bool sample::write()
{
    PtrList<scalarList> data;
    wordList headers;

    forAll(fields_, i)
    {
        this->appendData(fields_[i], data, headers);
    }

    if (Pstream::master() && (runTime_.value() > startTime_))
    {
        if (timeAverage_)
        {
            forAll(timeAveragedData_, i)
            {
                forAll(timeAveragedData_[i], j)
                {
                    timeAveragedData_[i][j]
                        += data[i][j] * runTime_.deltaTValue();
                }
            }

            if (runTime_.writeTime())
            {
                writeTime_ = runTime_.timeName();
            }

            const fileName path("postProcessing"/name_/writeTime_);
            mkDir(path);

            OFstream file(path/"averagedSample.txt");

            // Write header

            file<< "# x y z";

            forAll(headers, i)
            {
                file<< " " << headers[i];
            }

            file<< nl;

            // Write data

            const vectorList& points = interpPtr_->points();

            forAll(points, i)
            {
                file<< points[i].x() << " "
                    << points[i].y() << " "
                    << points[i].z();

                forAll(timeAveragedData_, j)
                {
                    const scalar avg
                        = timeAveragedData_[j][i]
                        /
                        (
                            runTime_.value()
                            - startTime_
                        );

                    file<< " " << avg;
                }

                file<< nl;
            }
        }
        else
        {
            const fileName path("postProcessing"/name_/runTime_.timeName());
            mkDir(path);

            OFstream file(path/"sample.txt");

            // Write header

            file<< "# x y z";

            forAll(headers, i)
            {
                file<< " " << headers[i];
            }

            file<< nl;

            // Write data

            const vectorList& points = interpPtr_->points();

            forAll(points, i)
            {
                file<< points[i].x() << " "
                    << points[i].y() << " "
                    << points[i].z();

                forAll(data, j)
                {
                    file<< " " << data[j][i];
                }

                file<< nl;
            }
        }
    }

    return true;
}

}

}

}

}
