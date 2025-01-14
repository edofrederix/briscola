#include "EigenSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SolverType, int Order>
EigenSolver<SolverType,Order>::EigenSolver(const dictionary& dict)
:
    EigenSolverBase(dict)
{}

template<class SolverType, int Order>
EigenSolver<SolverType,Order>::EigenSolver(const EigenSolver& s)
:
    EigenSolverBase(s)
{}

template<class SolverType, int Order>
EigenSolver<SolverType,Order>::~EigenSolver()
{}

template<class SolverType, int Order>
void EigenSolver<SolverType,Order>::prepare(const MatrixType& A)
{
    solverPtr_.reset(new SolverType(A));
    solverPtr_->compute(A);
}

template<class SolverType, int Order>
void EigenSolver<SolverType,Order>::solve(RhsType& x, const RhsType& b)
{
    if (solverPtr_.valid())
    {
        x = solverPtr_->solve(b);
    }
    else
    {
        FatalErrorInFunction
            << "Solver not prepared" << endl << abort(FatalError);
    }
}

}

}

}
