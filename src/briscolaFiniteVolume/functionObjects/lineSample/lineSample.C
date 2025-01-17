#include "lineSample.H"
#include "Time.H"
#include "addToRunTimeSelectionTable.H"
#include "fvMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

defineTypeNameAndDebug(lineSample, 0);

addToRunTimeSelectionTable
(
    functionObject,
    lineSample,
    dictionary
);

lineSample::lineSample
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
    endPoints_(dict.lookupOrDefault<Switch>("endPoints", false))
{
    init();
}

lineSample::~lineSample()
{}

vectorList lineSample::points()
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
