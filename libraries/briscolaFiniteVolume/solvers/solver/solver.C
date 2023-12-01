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
List<Type> solver<SType,Type,MeshType>::normFactors
(
    const linearSystem<SType,Type,MeshType>& Eqn,
    const label l
) const
{
    const meshLevel<SType,MeshType>& A = Eqn.A()[l];

    const meshLevel<Type,MeshType>& x = Eqn.x()[l];
    const meshLevel<Type,MeshType>& b = Eqn.b()[l];

    const meshLevel<Type,MeshType> Ay(rowSum(A)*gAverage(x));

    return max
    (
        gSum(cmptMag(Amul(A,x) - Ay) + cmptMag(b - Ay)),
        1e-20*pTraits<Type>::one
    );
}

template<class SType, class Type, class MeshType>
void solver<SType,Type,MeshType>::RBGS
(
    linearSystem<SType,Type,MeshType>& sys,
    List<Type>& xi,
    const label l,
    const label sweeps,
    const labelList& converged,
    const scalar omega
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const meshLevel<SType,MeshType>& A = sys.A()[l];
    const meshLevel<Type,MeshType>& b = sys.b()[l];

    bool singular = xi.size() > 0;

    if (!singular)
        xi = List<Type>(x.size(), Zero);

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        if (singular)
            forAll(xi, d)
                xi[d] = gAverage(x[d]);

        forAll(x, d)
        if (!converged[d])
        {
            // Even points

            meshDirection<Type,MeshType>& xd = x[d];

            const meshDirection<SType,MeshType>& Ad = A[d];
            const meshDirection<Type,MeshType>& bd = b[d];

            forAllCells(xd, i, j, k)
            {
                if ((i+j+k) % 2 == 0)
                    xd(i,j,k) +=
                        omega
                      * (bd(i,j,k) - Amul(Ad,xd,i,j,k) - xi[d])
                      / Ad(i,j,k).center();
            }
        }

        x.correctCommBoundaryConditions();

        forAll(x, d)
        if (!converged[d])
        {
            // Odd points

            meshDirection<Type,MeshType>& xd = x[d];

            const meshDirection<SType,MeshType>& Ad = A[d];
            const meshDirection<Type,MeshType>& bd = b[d];

            forAllCells(xd, i, j, k)
            {
                if ((i+j+k) % 2 == 1)
                    xd(i,j,k) +=
                        omega
                      * (bd(i,j,k) - Amul(Ad,xd,i,j,k) - xi[d])
                      / Ad(i,j,k).center();
            }
        }

        x.correctCommBoundaryConditions();
    }

    x.correctNonCommBoundaryConditions();

    if (!singular)
        xi.clear();
}

template<class SType, class Type, class MeshType>
void solver<SType,Type,MeshType>::LEXGS
(
    linearSystem<SType,Type,MeshType>& sys,
    List<Type>& xi,
    const label l,
    const label sweeps,
    const labelList& converged,
    const scalar omega
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const meshLevel<SType,MeshType>& A = sys.A()[l];
    const meshLevel<Type,MeshType>& b = sys.b()[l];

    bool singular = xi.size() > 0;

    if (!singular)
        xi = List<Type>(x.size(), Zero);

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        if (singular)
            forAll(xi, d)
                xi[d] = gAverage(x[d]);

        forAll(x, d)
        if (!converged[d])
        {
            meshDirection<Type,MeshType>& xd = x[d];

            const meshDirection<SType,MeshType>& Ad = A[d];
            const meshDirection<Type,MeshType>& bd = b[d];

            forAllCells(xd, i, j, k)
            {
                xd(i,j,k) +=
                    omega
                  * (bd(i,j,k) - Amul(Ad,xd,i,j,k) - xi[d])
                  / Ad(i,j,k).center();
            }
        }

        x.correctCommBoundaryConditions();
    }

    x.correctNonCommBoundaryConditions();

    if (!singular)
        xi.clear();
}

template<class SType, class Type, class MeshType>
void solver<SType,Type,MeshType>::JAC
(
    linearSystem<SType,Type,MeshType>& sys,
    List<Type>& xi,
    const label l,
    const label sweeps,
    const labelList& converged,
    const scalar omega
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const meshLevel<SType,MeshType>& A = sys.A()[l];
    const meshLevel<Type,MeshType>& b = sys.b()[l];

    bool singular = xi.size() > 0;

    if (!singular)
        xi = List<Type>(x.size(), Zero);

    meshLevel<Type,MeshType> y(x);

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        if (singular)
            forAll(xi, d)
                xi[d] = gAverage(x[d]);

        forAll(x, d)
        if (!converged[d])
        {
            meshDirection<Type,MeshType>& yd = y[d];

            const meshDirection<Type,MeshType>& xd = x[d];
            const meshDirection<SType,MeshType>& Ad = A[d];
            const meshDirection<Type,MeshType>& bd = b[d];

            forAllCells(xd, i, j, k)
            {
                yd(i,j,k) +=
                    omega
                  * (bd(i,j,k) - Amul(Ad,xd,i,j,k) - xi[d])
                  / Ad(i,j,k).center();
            }
        }

        x = y;

        x.correctCommBoundaryConditions();
    }

    x.correctNonCommBoundaryConditions();

    if (!singular)
        xi.clear();
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
