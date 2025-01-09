#ifndef KSPBCGSPETScSolver_H
#define KSPBCGSPETScSolver_H

#include "PETScKSPSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class KSPBCGSPETScSolver
:
    public PETScKSPSolver
{
public:

    // Runtime type information
    TypeName("KSPBCGS");

    // Construct from dictionary
    KSPBCGSPETScSolver(const dictionary& dict);

    // Copy constructor
    KSPBCGSPETScSolver(const KSPBCGSPETScSolver&);

    // Destructor
    virtual ~KSPBCGSPETScSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm);

    // Parallel solver bool
    virtual bool parallel() const { return true; }
};

}

}

}

#ifdef NoRepository
    #include "KSPBCGSPETScSolver.C"
#endif

#endif
