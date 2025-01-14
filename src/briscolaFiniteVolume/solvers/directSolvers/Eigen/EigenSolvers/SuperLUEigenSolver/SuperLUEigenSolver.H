#ifndef SuperLUEigenSolver_H
#define SuperLUEigenSolver_H

#include <Eigen/SuperLUSupport>
#include "EigenSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

typedef ::Eigen::SuperLU<EigenSolverBase::ColMajorMatrixType>
    SuperLUSolver;

class SuperLUEigenSolver
:
    public EigenSolver<SuperLUSolver, ::Eigen::ColMajor>
{
public:

    // Runtime type information
    TypeName("SuperLU");

    // Construct from dictionary
    SuperLUEigenSolver(const dictionary& dict);

    // Copy constructor
    SuperLUEigenSolver(const SuperLUEigenSolver&);

    // Destructor
    virtual ~SuperLUEigenSolver();
};

}

}

}

#ifdef NoRepository
    #include "SuperLUEigenSolver.C"
#endif

#endif
