#include "solver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
solver<SType,Type,MeshType>::solver
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    immersedBoundary<Type,MeshType>* IB
)
:
    dict_(dict),
    fvMsh_(fvMsh),
    relTol_(readScalar(dict.lookup("relTol"))),
    tolerance_(readScalar(dict.lookup("tolerance"))),
    minIter_(dict.lookupOrDefault<label>("minIter", 0)),
    maxIter_(dict.lookupOrDefault<label>("maxIter", 999)),
    IB_(IB)
{}

template<class SType, class Type, class MeshType>
solver<SType,Type,MeshType>::~solver()
{}

template<class SType, class Type, class MeshType>
autoPtr<solver<SType,Type,MeshType>> solver<SType,Type,MeshType>::New
(
    const word solverName,
    const fvMesh& fvMsh,
    immersedBoundary<Type,MeshType>* IB
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

    return autoPtr<solver<SType,Type,MeshType>>(cstrIter()(dict, fvMsh, IB));
}

template<class SType, class Type, class MeshType>
autoPtr<solver<SType,Type,MeshType>> solver<SType,Type,MeshType>::New
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    immersedBoundary<Type,MeshType>* IB
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

    return autoPtr<solver<SType,Type,MeshType>>(cstrIter()(dict, fvMsh, IB));
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
    const label l,
    const label sweeps,
    const labelList& converged,
    const scalar omega,
    immersedBoundary<Type, MeshType>* IB
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const meshLevel<SType,MeshType>& A = sys.A()[l];
    const meshLevel<Type,MeshType>& b = sys.b()[l];

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
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
                      * (bd(i,j,k) - Amul(Ad,xd,i,j,k))
                      / Ad(i,j,k).center();
            }
        }

        x.correctCommBoundaryConditions();
        correctImmersedBoundaryConditions(IB,sys,l);

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
                      * (bd(i,j,k) - Amul(Ad,xd,i,j,k))
                      / Ad(i,j,k).center();
            }
        }

        x.correctCommBoundaryConditions();
        correctImmersedBoundaryConditions(IB,sys,l);
    }

    x.correctNonCommBoundaryConditions();
}

template<class SType, class Type, class MeshType>
void solver<SType,Type,MeshType>::LEXGS
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label sweeps,
    const labelList& converged,
    const scalar omega,
    immersedBoundary<Type, MeshType>* IB
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const meshLevel<SType,MeshType>& A = sys.A()[l];
    const meshLevel<Type,MeshType>& b = sys.b()[l];

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
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
                  * (bd(i,j,k) - Amul(Ad,xd,i,j,k))
                  / Ad(i,j,k).center();
            }
        }

        x.correctCommBoundaryConditions();
        correctImmersedBoundaryConditions(IB,sys,l);
    }

    x.correctNonCommBoundaryConditions();
}

template<class SType, class Type, class MeshType>
void solver<SType,Type,MeshType>::JAC
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label sweeps,
    const labelList& converged,
    const scalar omega,
    immersedBoundary<Type, MeshType>* IB
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const meshLevel<SType,MeshType>& A = sys.A()[l];
    const meshLevel<Type,MeshType>& b = sys.b()[l];

    meshLevel<Type,MeshType> y(x);

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
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
                  * (bd(i,j,k) - Amul(Ad,xd,i,j,k))
                  / Ad(i,j,k).center();
            }
        }

        x = y;

        x.correctCommBoundaryConditions();
        correctImmersedBoundaryConditions(IB,sys,l);
    }

    x.correctNonCommBoundaryConditions();
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
void solver<SType,Type,MeshType>::correctImmersedBoundaryConditions
(
    immersedBoundary<Type, MeshType>* IB,
    linearSystem<SType,Type,MeshType>& sys,
    const label l
)
{
    if (IB != nullptr)
    {
        meshLevel<Type,MeshType>& x = sys.x()[l];

        // Penalization
        // forAllDirections(x,d,i,j,k)
        // {
        //     if (Foam::mag(IB->mask()(l,d,i,j,k)-1.0) < 0.01)
        //     {
        //         x(d,i,j,k) = Zero;
        //     }
        // }

        // Fadlun method
        // forAllDirections(x,d,i,j,k)
        // {
        //     const labelVector ijk(i,j,k);
        //     if (Foam::mag(IB->wallAdjMask()(l,d,i,j,k)-1.0) < 0.01)
        //     {
        //         scalar ximax = 0;
        //         for (int dir = 0; dir < 6; dir++)
        //         {
        //             const label oppositeDir =
        //                 faceNumber(-faceOffsets[dir]);

        //             if (IB->wallDist()(l,d,i,j,k)[dir] > ximax)
        //             {
        //                 ximax = IB->wallDist()(l,d,i,j,k)[dir];
        //                 const scalar xic = 1.0 - IB->wallDist()(l,d,i,j,k)[dir];
        //                 const scalar xinb = 1.0 + xic;
        //                 const scalar w = xic/xinb;

        //                 x(d,i,j,k) = w*x[d](ijk+faceOffsets[oppositeDir]);
        //             }
        //         }
        //     }
        // }

        // Deen method
        forAllDirections(x,d,i,j,k)
        {
            const labelVector ijk(i,j,k);
            if (Foam::mag(IB->wallAdjMask()(l,d,i,j,k)-1.0) < 0.01)
            {
                for (int dir = 0; dir < 6; dir++)
                {
                    const label oppositeDir =
                        faceNumber(-faceOffsets[dir]);

                    if (IB->wallDist()(l,d,i,j,k)[dir] > 0)
                    {
                        scalar xiStabilityFactor = 0;

                        const scalar xi
                            = IB->wallDist()(l,d,i,j,k)[dir];

                        const scalar xi2
                            = IB->neighborDist()(l,d,i,j,k)[oppositeDir];

                        scalar w0 = 2.0 /
                            (
                                (1.0 - xiStabilityFactor)
                                * (2.0 - xiStabilityFactor)
                            );
                        scalar w1 = 2.0 - (2.0 - xi) * w0;
                        scalar w2 = -1.0 + (1.0 - xi) * w0;

                        if (xi < xiStabilityFactor)
                        {
                            w1 = xi*xi2/((1.0-xi)*(1.0-xi2));
                            w2 = xi/((xi2-xi)*(xi2-1.0));
                        }

                        x[d](ijk+faceOffsets[dir])
                            = w1*x[d](ijk) + w2*x[d](ijk-faceOffsets[dir]);
                    }
                }
            }
        }
    }
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
