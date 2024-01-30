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

bool EigenSolver::isPositiveDefinite(const matrixType& A) const
{
    if (!isSymmetric(A))
        return false;

    // The check below is sufficient in practice

    for (int k = 0; k < A.outerSize(); ++k)
    {
        scalar C = Zero;

        for (matrixType::InnerIterator it(A,k); it; ++it)
        {
            const int i = it.row();
            const int j = it.col();

            if (i == j)
            {
                C += Foam::mag(A.coeff(i,j));
            }
            else if (i < j)
            {
                C -= 2.0*Foam::mag(A.coeff(i,j));
            }
        }

        if (C < -1e-12)
            return false;
    }

    return true;
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
