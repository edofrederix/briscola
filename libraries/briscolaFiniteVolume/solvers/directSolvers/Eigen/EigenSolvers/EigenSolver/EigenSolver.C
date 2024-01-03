#include "EigenSolver.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(EigenSolver, 0);
defineRunTimeSelectionTable(EigenSolver, dictionary);

bool EigenSolver::isSymmetric(const matrixType& A) const
{
    bool symm = true;

    // Only check non-zero elements

    for (int k = 0; k < A.outerSize(); ++k)
    for (matrixType::InnerIterator it(A,k); it; ++it)
    {
        const int i = it.row();
        const int j = it.col();

        if(i < j && Foam::mag(A.coeff(i,j) - A.coeff(j,i)) > 1e-8)
        {
            symm = false;
            goto done;
        }
    }

    done:;

    return symm;
}

EigenSolver::EigenSolver()
{}

EigenSolver::EigenSolver(const EigenSolver& s)
{}

EigenSolver::~EigenSolver()
{}

autoPtr<EigenSolver> EigenSolver::New(const word solverType)
{
    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown Eigen solver type " << solverType
            << ". Valid Eigen solvers are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<EigenSolver>(cstrIter()());
}

}

}

}
