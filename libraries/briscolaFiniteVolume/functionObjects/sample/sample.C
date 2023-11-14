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
    name_(name)
{}

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

    if (Pstream::master())
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

    return true;
}

}

}

}

}
