#ifndef functionObjects_timeAverage_H
#define functionObjects_timeAverage_H

#include "briscolaFunctionObject.H"
#include "fvMesh.H"
#include "meshField.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace functionObjects
{

class timeAverage
:
    public briscolaFunctionObject
{

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Fields to be time averaged
    wordList fields_;

    // Time averaged colocated scalar fields
    PtrList<meshField<scalar,colocated>> avgColScalarFields_;

    // Time averaged colocated vector fields
    PtrList<meshField<vector,colocated>> avgColVectorFields_;

    // Time averaged colocated tensor fields
    PtrList<meshField<tensor,colocated>> avgColTensorFields_;

    // Time averaged colocated symm tensor fields
    PtrList<meshField<symmTensor,colocated>> avgColSymmTensorFields_;

    // Time averaged staggered scalar fields
    PtrList<meshField<scalar,staggered>> avgStagScalarFields_;

    // Time averaged staggered vector fields
    PtrList<meshField<vector,staggered>> avgStagVectorFields_;

    // Averaging start time
    scalar startTime_;

    // Average in intervals
    Switch intervalAveraging_;

    // Start time of averaging interval
    scalar intervalStart_;

    // Switch to reset averages for interval averaging
    Switch reset_;

    // Initialize
    void init();

public:

    // Runtime type information
    TypeName("timeAverage");

    // Construct from time and dictionary
    timeAverage
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict
    );

    // Disallow default bitwise copy construction
    timeAverage(const timeAverage&) = delete;

    // Destructor
    virtual ~timeAverage();

    // Code executed at the end of each time step loop
    virtual bool execute();

    // Code executed at the end of each time step loop, after execute()
    virtual bool write();

    // Code executed at the end of the simulation
    virtual bool end();

    // Disallow default bitwise assignment
    void operator=(const timeAverage&) = delete;
};

}

}

}

}

#endif
