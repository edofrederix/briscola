#ifndef functionObjects_lineSample_H
#define functionObjects_lineSample_H

#include "sample.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class lineSample
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

public:

    // Runtime type information
    TypeName("lineSample");

    // Construct from time and dictionary
    lineSample
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    lineSample(const lineSample&) = delete;

    // Destructor
    virtual ~lineSample();

    // Derived functions should define the points
    virtual vectorList points();

    // Disallow default bitwise assignment
    void operator=(const lineSample&) = delete;
};

}

}

}

}

#endif
