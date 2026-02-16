#include "averagedLineSample.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"
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

defineTypeNameAndDebug(averagedLineSample, 0);

addToRunTimeSelectionTable
(
    functionObject,
    averagedLineSample,
    dictionary
);

averagedLineSample::averagedLineSample
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    sample(name, runTime, dict),
    start_(dict.lookup("start")),
    end_(dict.lookup("end")),
    N_(readLabel(dict.lookup("N"))),
    endPoints_(dict.lookupOrDefault<Switch>("endPoints", false)),
    averagingDirections_(dict.lookup("averagingDirections"))
{
    for (int dir = 0; dir < 3; dir++)
    {
        if
        (
               (averagingDirections_[dir] != 1)
            && (averagingDirections_[dir] != 0)
        )
        {
            FatalErrorInFunction
                << "Invalid averaging direction vector." << endl
                << abort(FatalError);
        }
    }

    /* Check that averaging directions are perpendicular to line direction */

    init();

    /* Check that time-averaging is on */

    const objectRegistry& db = this->fvMsh_.db();

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

    averagedData_.setSize(size);

    forAll(averagedData_, i)
    {
        averagedData_.set
        (
            i,
            new scalarList
            (
                N_,
                Zero
            )
        );
    }
}

averagedLineSample::~averagedLineSample()
{}

bool averagedLineSample::write()
{
    return true;
}

bool averagedLineSample::end()
{
    // Mesh
    const mesh& msh = this->fvMsh_.msh();

    // Cell sizes
    const FastPtrList<PartialList<scalar>>& cellSizes =
        msh.cast<rectilinearMesh>().globalCellSizes();

    vectorList points = this->points();

    labelVector averagingN(1,1,1);

    forAll(averagedData_, i)
    {
        forAll(averagedData_[i], j)
        {
            averagedData_[i][j] = Zero;
        }
    }

    wordList headersList;

    // set line to bounding box edge in directions that we are averaging
    if (averagingDirections_.x())
    {
        forAll(points, i)
        {
            points[i].x() = msh[0].boundingBox().left();
        }

        averagingN.x() = msh.cast<rectilinearMesh>().N().x();
    }
    if (averagingDirections_.y())
    {
        forAll(points, i)
        {
            points[i].y() = msh[0].boundingBox().bottom();
        }

        averagingN.y() = msh.cast<rectilinearMesh>().N().y();
    }
    if (averagingDirections_.z())
    {
        forAll(points, i)
        {
            points[i].z() = msh[0].boundingBox().aft();
        }

        averagingN.z() = msh.cast<rectilinearMesh>().N().z();
    }

    // Main spatial averaging loops
    for (int i = 0; i < averagingN.x(); i++)
    {
        if (averagingDirections_.y())
        {
            forAll(points, j)
            {
                points[j].y() = msh[0].boundingBox().bottom();
            }
        }

        for (int j = 0; j < averagingN.y(); j++)
        {
            if (averagingDirections_.z())
            {
                forAll(points, l)
                {
                    points[l].z() = msh[0].boundingBox().aft();
                }
            }

            for (int k = 0; k < averagingN.z(); k++)
            {
                this->interpPtr_.reset
                (
                    pointInterpolator<colocated>::New
                    (
                        this->fvMsh_,
                        points,
                        this->dict_.lookupOrDefault<word>("interpolator", "linear")
                    ).ptr()
                );

                FastPtrList<scalarList> data;
                wordList headers;

                forAll(fields_, l)
                {
                    this->appendData(fields_[l], data, headers);
                }

                headersList = headers;

                forAll(averagedData_, l)
                {
                    forAll(averagedData_[l], m)
                    {
                        averagedData_[l][m] += data[l][m];
                    }
                }

                if (averagingDirections_.z())
                {
                    forAll(points, l)
                    {
                        points[l].z() += cellSizes[2][k];
                    }
                }
            }

            if (averagingDirections_.y())
            {
                forAll(points, l)
                {
                    points[l].y() += cellSizes[1][j];
                }
            }
        }

        if (averagingDirections_.x())
        {
            forAll(points, l)
            {
                points[l].x() += cellSizes[0][i];
            }
        }
    }

    forAll(averagedData_, i)
    {
        forAll(averagedData_[i], j)
        {
            averagedData_[i][j] /= cmptProduct(averagingN);
        }
    }

    if (Pstream::master() && (runTime_.value() > startTime_))
    {
        const fileName path("postProcessing"/name_);
        mkDir(path);

        OFstream file(path/"averagedSample.txt");

        // Write header

        file<< "# x y z";

        forAll(headersList, i)
        {
            file<< " " << headersList[i];
        }

        file<< nl;

        // Write data

        points = this->points();

        forAll(points, i)
        {
            file<< points[i].x() << " "
                << points[i].y() << " "
                << points[i].z();

            forAll(averagedData_, j)
            {
                const scalar avg
                    = averagedData_[j][i]
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

    return true;
}

vectorList averagedLineSample::points()
{
    vectorList points(N_);

    forAll(points, i)
    {
        if (endPoints_)
        {
            points[i] = start_ + (end_ - start_)*i/(N_-1);
        }
        else
        {
            points[i] = start_ + (end_ - start_)*(i+0.5)/N_;
        }
    }

    return points;
}

}

}

}

}
