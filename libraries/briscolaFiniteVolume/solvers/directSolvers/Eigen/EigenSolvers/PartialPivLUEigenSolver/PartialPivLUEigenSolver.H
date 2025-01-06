#ifndef PartialPivLUEigenSolver_H
#define PartialPivLUEigenSolver_H

#include <Eigen/Dense>
#include "EigenSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

typedef ::Eigen::PartialPivLU<EigenSolverBase::RhsType>
    PartialPivLUSolver;

class PartialPivLUEigenSolver
:
    public EigenSolver<PartialPivLUSolver, ::Eigen::ColMajor>
{
public:

    // Runtime type information
    TypeName("PartialPivLU");

    // Construct from dictionary
    PartialPivLUEigenSolver(const dictionary& dict);

    // Copy constructor
    PartialPivLUEigenSolver(const PartialPivLUEigenSolver&);

    // Destructor
    virtual ~PartialPivLUEigenSolver();

    // First guess is not needed
    virtual bool needGuess() const { return false; }
};

}

}

}

#ifdef NoRepository
    #include "PartialPivLUEigenSolver.C"
#endif

#endif
