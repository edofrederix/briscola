#ifndef SparseLUEigenSolver_H
#define SparseLUEigenSolver_H

#include "EigenSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

typedef ::Eigen::SparseLU<EigenSolverBase::ColMajorMatrixType>
    SparseLUSolver;

class SparseLUEigenSolver
:
    public EigenSolver<SparseLUSolver, ::Eigen::ColMajor>
{
public:

    // Runtime type information
    TypeName("SparseLU");

    // Construct from dictionary
    SparseLUEigenSolver(const dictionary& dict);

    // Copy constructor
    SparseLUEigenSolver(const SparseLUEigenSolver&);

    // Destructor
    virtual ~SparseLUEigenSolver();

    // First guess is not needed
    virtual bool needGuess() const { return false; }
};

}

}

}

#ifdef NoRepository
    #include "SparseLUEigenSolver.C"
#endif

#endif
