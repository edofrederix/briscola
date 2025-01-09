#ifndef KSPGMRESPETScSolver_H
#define KSPGMRESPETScSolver_H

#include "PETScKSPSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

class KSPGMRESPETScSolver
:
    public PETScKSPSolver
{
public:

    // Runtime type information
    TypeName("KSPGMRES");

    // Construct from dictionary
    KSPGMRESPETScSolver(const dictionary& dict);

    // Copy constructor
    KSPGMRESPETScSolver(const KSPGMRESPETScSolver&);

    // Destructor
    virtual ~KSPGMRESPETScSolver();

    // Prepare from a sparse matrix
    virtual void prepare(const Mat& mat, const MPI_Comm& comm);

    // Parallel solver bool
    virtual bool parallel() const { return true; }
};

}

}

}

#ifdef NoRepository
    #include "KSPGMRESPETScSolver.C"
#endif

#endif
