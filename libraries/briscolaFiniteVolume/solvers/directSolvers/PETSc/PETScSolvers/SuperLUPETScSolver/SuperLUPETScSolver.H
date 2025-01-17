#ifndef SuperLUPETScSolver_H
#define SuperLUPETScSolver_H

#include "PETScDirectSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class SuperLUPETScSolver
:
    public PETScDirectSolver
{
public:

    // Runtime type information
    TypeName("SuperLU");

    // Construct from dictionary
    SuperLUPETScSolver(const dictionary& dict);

    // Copy constructor
    SuperLUPETScSolver(const SuperLUPETScSolver&);

    // Destructor
    virtual ~SuperLUPETScSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm);

    // Parallel solver bool
    virtual bool parallel() const { return false; }
};

}

}

}

#ifdef NoRepository
    #include "SuperLUPETScSolver.C"
#endif

#endif
