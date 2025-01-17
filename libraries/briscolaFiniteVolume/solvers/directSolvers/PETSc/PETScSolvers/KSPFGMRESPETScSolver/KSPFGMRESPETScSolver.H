#ifndef KSPFGMRESPETScSolver_H
#define KSPFGMRESPETScSolver_H

#include "PETScKSPSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class KSPFGMRESPETScSolver
:
    public PETScKSPSolver
{
public:

    // Runtime type information
    TypeName("KSPFGMRES");

    // Construct from dictionary
    KSPFGMRESPETScSolver(const dictionary& dict);

    // Copy constructor
    KSPFGMRESPETScSolver(const KSPFGMRESPETScSolver&);

    // Destructor
    virtual ~KSPFGMRESPETScSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm);

    // Parallel solver bool
    virtual bool parallel() const { return true; }
};

}

}

}

#ifdef NoRepository
    #include "KSPFGMRESPETScSolver.C"
#endif

#endif
