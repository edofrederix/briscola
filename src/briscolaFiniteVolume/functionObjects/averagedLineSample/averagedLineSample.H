#ifndef functionObjects_averagedLineSample_H
#define functionObjects_averagedLineSample_H

#include "sample.H"
#include "List.H"
#include "rectilinearMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class averagedLineSample
:
    public sample
{
protected:

    // Start point
    vector start_;

    // End point
    vector end_;

    // Number of sample points
    label N_;

    // Include end points. Defaults to false.
    bool endPoints_;

    // Directions in which spatial averaging is applied
    labelVector averagingDirections_;

    // Spatially averaged data
    PtrList<scalarList> averagedData_;

public:

    // Runtime type information
    TypeName("averagedLineSample");

    // Construct from time and dictionary
    averagedLineSample
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    averagedLineSample(const averagedLineSample&) = delete;

    // Destructor
    virtual ~averagedLineSample();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Code executed at the end of the simulation
    virtual bool end();

    // Derived functions should define the points
    virtual vectorList points();

    // Disallow default bitwise assignment
    void operator=(const averagedLineSample&) = delete;
};

}

}

}

}

#endif
