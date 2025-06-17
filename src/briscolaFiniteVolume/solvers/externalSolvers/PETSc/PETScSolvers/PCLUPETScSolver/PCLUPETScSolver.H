#ifndef PCLUPETScSolver_H
#define PCLUPETScSolver_H

#include "PETScExternalSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class PCLUPETScSolver
:
    public PETScExternalSolver
{
public:

    // Runtime type information
    TypeName("PCLU");

    // Construct from dictionary
    PCLUPETScSolver(const dictionary& dict);

    // Copy constructor
    PCLUPETScSolver(const PCLUPETScSolver&);

    // Destructor
    virtual ~PCLUPETScSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm);

    // Parallel solver bool
    virtual bool parallel() const { return false; }
};

}

}

}

#ifdef NoRepository
    #include "PCLUPETScSolver.C"
#endif

#endif
