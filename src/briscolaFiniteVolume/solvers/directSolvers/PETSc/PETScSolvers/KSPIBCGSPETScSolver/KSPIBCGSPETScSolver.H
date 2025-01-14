#ifndef KSPIBCGSPETScSolver_H
#define KSPIBCGSPETScSolver_H

#include "PETScKSPSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class KSPIBCGSPETScSolver
:
    public PETScKSPSolver
{
public:

    // Runtime type information
    TypeName("KSPIBCGS");

    // Construct from dictionary
    KSPIBCGSPETScSolver(const dictionary& dict);

    // Copy constructor
    KSPIBCGSPETScSolver(const KSPIBCGSPETScSolver&);

    // Destructor
    virtual ~KSPIBCGSPETScSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm);

    // Parallel solver bool
    virtual bool parallel() const { return true; }
};

}

}

}

#ifdef NoRepository
    #include "KSPIBCGSPETScSolver.C"
#endif

#endif
