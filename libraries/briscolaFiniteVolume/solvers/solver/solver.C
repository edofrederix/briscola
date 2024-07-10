#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
solver<SType,Type,MeshType>::solver(const dictionary& dict, const fvMesh& fvMsh)
:
    dict_(dict),
    fvMsh_(fvMsh),
    relTol_(readScalar(dict.lookup("relTol"))),
    tolerance_(readScalar(dict.lookup("tolerance"))),
    minIter_(dict.lookupOrDefault<label>("minIter", 0)),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 999))
{}

template<class SType, class Type, class MeshType>
solver<SType,Type,MeshType>::~solver()
{}

template<class SType, class Type, class MeshType>
autoPtr<solver<SType,Type,MeshType>> solver<SType,Type,MeshType>::New
(
    const word solverName,
    const fvMesh& fvMsh
)
{
    const dictionary dict
    (
        fvMsh.solverDict().subDict(solverName)
    );

    const word solverType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown solver "
            << solverType << nl << nl
            << "Valid solvers are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<solver<SType,Type,MeshType>>(cstrIter()(dict, fvMsh));
}

template<class SType, class Type, class MeshType>
autoPtr<solver<SType,Type,MeshType>> solver<SType,Type,MeshType>::New
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
{
    const word solverType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(solverType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown solver "
            << solverType << nl << nl
            << "Valid solvers are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<solver<SType,Type,MeshType>>(cstrIter()(dict, fvMsh));
}

template<class SType, class Type, class MeshType>
void solver<SType,Type,MeshType>::solve
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const bool constMatrix
)
{
    this->solve
    (
        const_cast<tmp<linearSystem<SType,Type,MeshType>>&>(tSys).ref(),
        constMatrix
    );

    tSys.clear();
}

template<class SType, class Type, class MeshType>
List<Type> solver<SType,Type,MeshType>::normFactors
(
    linearSystem<SType,Type,MeshType>& sys,
    const meshField<Type,MeshType>& res,
    const label l
) const
{
    List<Type> y(gAverage(sys.x()[l]));
    List<Type> f(MeshType::numberOfDirections, Zero);

    const List<bool> diagonal(sys.diagonal());

    forAll(f, d)
    {
        const meshDirection<SType,MeshType>& A = sys.A()[l][d];
        const meshDirection<Type,MeshType>& b = sys.b()[l][d];
        const meshDirection<Type,MeshType>& r = res[l][d];

        forAllCells(A, i, j, k)
        {
            Type Ay =
                diagonal[d]
              ? A(i,j,k).center()*y[d]
              : rowSum(A,i,j,k)*y[d];

            f[d] +=
                Foam::cmptMag(r(i,j,k) - b(i,j,k) + Ay)
              + Foam::cmptMag(b(i,j,k) - Ay);
        }
    }

    reduce(f, sumOp<List<Type>>());

    return max(f, 1e-20*pTraits<Type>::one);
}

template<class SType, class Type, class MeshType>
void solver<SType,Type,MeshType>::printSolverStats
(
    const word solverName,
    const word fieldName,
    const List<Type> initialResiduals,
    const List<Type> finalResiduals,
    const label nIter
)
{
    Info<< solverName
        << ": Solving for " << MeshType::typeName << " " << fieldName
        << ", initial residual = " << max(cmptMax(initialResiduals))
        << ", final residual = " << max(cmptMax(finalResiduals))
        << ", nIter = " << nIter << endl;
}

template<class SType, class Type, class MeshType>
labelList solver<SType,Type,MeshType>::checkConvergence
(
    const List<Type>& currentResiduals,
    const List<Type>& initialResiduals,
    const scalar relTol,
    const scalar absTol
) const
{
    labelList conv(MeshType::numberOfDirections);

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        conv[d] =
            (
                (Foam::cmptMax(initialResiduals[d]) > 0.0)
             && (Foam::cmptMax(currentResiduals[d]) > absTol)
             && (
                    Foam::cmptMax
                    (
                        Foam::cmptDivide
                        (
                            currentResiduals[d],
                            initialResiduals[d] + pTraits<Type>::one*1e-20
                        )
                    )
                  > relTol
                )
            ) ? 0 : 1;
    }

    return conv;
}

}

}

}
