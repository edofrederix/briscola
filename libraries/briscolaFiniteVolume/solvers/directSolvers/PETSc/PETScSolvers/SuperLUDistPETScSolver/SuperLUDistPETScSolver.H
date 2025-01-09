#ifndef SuperLUDistPETScSolver_H
#define SuperLUDistPETScSolver_H

#include "PETScDirectSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class SuperLUDistPETScSolver
:
    public PETScDirectSolver
{
public:

    // Runtime type information
    TypeName("SuperLUDist");

    // Construct from dictionary
    SuperLUDistPETScSolver(const dictionary& dict);

    // Copy constructor
    SuperLUDistPETScSolver(const SuperLUDistPETScSolver&);

    // Destructor
    virtual ~SuperLUDistPETScSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm);

    // Parallel solver bool
    virtual bool parallel() const { return true; }
};

}

}

}

#ifdef NoRepository
    #include "SuperLUDistPETScSolver.C"
#endif

#endif
