#ifndef functionObjects_sample_H
#define functionObjects_sample_H

#include "briscolaFunctionObject.H"
#include "pointInterpolator.H"
#include "fvMesh.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class sample
:
    public briscolaFunctionObject
{
protected:

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Point interpolator
    autoPtr<pointInterpolator<colocated>> interpPtr_;

    // Fields to be sampled
    wordList fields_;

    // Name of the function object
    word name_;

    // Time average output data
    bool timeAverage_;

    // Time averaged data
    PtrList<scalarList> timeAveragedData_;

    // Sampling start time
    scalar startTime_;

    // Name of write-time folder
    word writeTime_;

    // Initialize
    void init();

    // Append to data array
    void appendData(const word, PtrList<scalarList>&, wordList&);

    // Templated append function

    template<class Type>
    void appendData
    (
        const meshField<Type,colocated>&,
        PtrList<scalarList>&,
        wordList&
    );

    template<class Type>
    void appendScalarData
    (
        const meshField<Type,colocated>&,
        PtrList<scalarList>&,
        wordList&
    );

    template<class Type>
    void appendArrayData
    (
        const meshField<Type,colocated>&,
        PtrList<scalarList>&,
        wordList&
    );

    template<class Type>
    void appendArrayArrayData
    (
        const meshField<Type,colocated>&,
        PtrList<scalarList>&,
        wordList&
    );

public:

    // Construct from time and dictionary
    sample
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    sample(const sample&) = delete;

    // Destructor
    virtual ~sample();

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Derived functions should define the points
    virtual vectorList points() = 0;

    // Disallow default bitwise assignment
    void operator=(const sample&) = delete;
};

}

}

}

}

#endif
