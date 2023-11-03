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
                if (((i+j+k) % 2 == 0) && (IB == nullptr ? true : (IB->ghostMask()(l,d,i,j,k) != 1)))
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
                if (((i+j+k) % 2 == 1) && (IB == nullptr ? true : (IB->ghostMask()(l,d,i,j,k) != 1)))
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
                if (IB == nullptr ? true : (IB->ghostMask()(l,d,i,j,k) != 1))
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
                if (IB == nullptr ? true : (IB->ghostMask()(l,d,i,j,k) != 1))
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
        scalar omega = 0.8;

        meshLevel<Type,MeshType>& x = sys.x()[l];

        // Penalization
        // forAllDirections(x,d,i,j,k)
        // {
        //     if (IB->mask()(l,d,i,j,k) == 1)
        //     {
        //         x(d,i,j,k) = Zero;
        //     }
        // }

        // Fadlun method
        // forAllDirections(x,d,i,j,k)
        // {
        //     const labelVector ijk(i,j,k);
        //     if (IB->wallAdjMask()(l,d,i,j,k) == 1)
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

        // Deen/Vreman method
        if (IB->type() == "Vreman")
        {
        forAllDirections(x,d,i,j,k)
        {
            const labelVector ijk(i,j,k);
            if (IB->wallAdjMask()(l,d,i,j,k) == 1)
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

                        x[d](ijk+faceOffsets[dir]) = (1.0 - omega) * x(d,i,j,k)
                            + omega * (w1*x[d](ijk) + w2*x[d](ijk-faceOffsets[dir]));


                    }
                }
            }
        }
        }

        // Mittal

        if (IB->type() ==  "Mittal")
        {
        // Cell centers
        const fvMesh& fvMsh = sys.fvMsh();
        const meshField<vector,MeshType>& CC =
            fvMsh.metrics<MeshType>().cellCenters();
        const mesh& msh = fvMsh.msh();

        scalar tol = 1e-5;
        // Info << msh[l].boundingBox().fore() << endl;

        forAllDirections(x,d,i,j,k)
        {
            if (IB->ghostMask()(l,d,i,j,k) == 1)
            {
                // mirror point - make this a field in IB
                vector mp = IB->mirrorPoint(CC(l,d,i,j,k));

                // Fix situations where the mirror point is just outside of the mesh
                // bounding box due to rounding errors
                if
                (
                    (mp.x() <= msh[l].boundingBox().left() + tol)
                    && (mp.x() >= msh[l].boundingBox().left() - tol)
                )
                {
                    mp.x() += tol;
                }
                if
                (
                    (mp.x() >= msh[l].boundingBox().right() - tol)
                    && (mp.x() <= msh[l].boundingBox().right() + tol)
                )
                {
                    mp.x() -= tol;
                }
                if
                (
                    (mp.y() <= msh[l].boundingBox().bottom() + tol)
                    && (mp.y() >= msh[l].boundingBox().bottom() - tol)
                )
                {
                    mp.y() += tol;
                }
                if
                (
                    (mp.y() >= msh[l].boundingBox().top() - tol)
                    && (mp.y() <= msh[l].boundingBox().top() + tol)
                )
                {
                    mp.y() -= tol;
                }
                if
                (
                    (mp.z() <= msh[l].boundingBox().aft() + tol)
                    && (mp.z() >= msh[l].boundingBox().aft() - tol)
                )
                {
                    mp.z() += tol;
                }
                if
                (
                    (mp.z() >= msh[l].boundingBox().fore() - tol)
                    && (mp.z() <= msh[l].boundingBox().fore() + tol)
                )
                {
                    mp.z() -= tol;
                }

                // Colocated cell index of mp
                labelVector mpIndex = msh.findCell(mp, l);

                // Local coordinates of mp in colocated cell
                vector mpLocalCoords = msh[l].points().cellCoordinates(mp, mpIndex, true);

                if (mpLocalCoords == -vector::one)
                {
                    FatalError
                        << "Interpolation error at level " << l
                        << " and direction " << d
                        << ". Mirror point: " << mp
                        << " and colocated cell index: "
                        << mpIndex
                        << endl;
                    FatalError.exit();
                }

                // Index of staggered left-bottom-aft cell w.r.t. mp
                labelVector mpLBA = mpIndex;
                if
                (
                    (d != 0)
                    && (mpLocalCoords.x() < 0.5)
                )
                {
                    mpLBA.x() -= 1;
                }
                if
                (
                    (d != 1)
                    && (mpLocalCoords.y() < 0.5)
                )
                {
                    mpLBA.y() -= 1;
                }
                if
                (
                    (d != 2)
                    && (mpLocalCoords.z() < 0.5)
                )
                {
                    mpLBA.z() -= 1;
                }

                // Interpolation box
                vertexVector interpPoints
                (
                    CC[l][d](mpLBA),
                    CC[l][d](mpLBA+unitX),
                    CC[l][d](mpLBA+unitY),
                    CC[l][d](mpLBA+unitXY),
                    CC[l][d](mpLBA+unitZ),
                    CC[l][d](mpLBA+unitXZ),
                    CC[l][d](mpLBA+unitYZ),
                    CC[l][d](mpLBA+unitXYZ)
                );

                // Interpolation weights
                const vector v(interpolationWeights(mp,interpPoints,true));

                vertexScalar weights
                (
                    (1-v.x())*(1-v.y())*(1-v.z()),
                    (  v.x())*(1-v.y())*(1-v.z()),
                    (1-v.x())*(  v.y())*(1-v.z()),
                    (  v.x())*(  v.y())*(1-v.z()),
                    (1-v.x())*(1-v.y())*(  v.z()),
                    (  v.x())*(1-v.y())*(  v.z()),
                    (1-v.x())*(  v.y())*(  v.z()),
                    (  v.x())*(  v.y())*(  v.z())
                );

                Type mpValue = Zero;

                if (v != -vector::one)
                {
                    mpValue =
                          weights.lba()*x[d](mpLBA)
                        + weights.rba()*x[d](mpLBA+unitX)
                        + weights.lta()*x[d](mpLBA+unitY)
                        + weights.rta()*x[d](mpLBA+unitX+unitY)
                        + weights.lbf()*x[d](mpLBA+unitZ)
                        + weights.rbf()*x[d](mpLBA+unitX+unitZ)
                        + weights.ltf()*x[d](mpLBA+unitY+unitZ)
                        + weights.rtf()*x[d](mpLBA+unitX+unitY+unitZ);
                }
                else
                {
                    FatalError
                        << "Interpolation error at level " << l
                        << " and direction " << d
                        << ". Mirror point: " << mp << nl
                        << "and interpolation points: "
                        << interpPoints << nl
                        << "Local colocated coordinates: "
                        << mpLocalCoords << nl
                        << "LBA cell index: " << mpLBA << nl
                        << "Colocated cell index" << mpIndex
                        << endl;
                    FatalError.exit();
                }

                x(d,i,j,k) = (1.0 - omega) * x(d,i,j,k) - omega * mpValue;
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
