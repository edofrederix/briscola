#include "APLU.H"

#include "SquareMatrix.H"
#include "LUscalarMatrix.H"
#include "meshDirectionStencilFunctions.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
APLU<SType,Type,MeshType>::APLU
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
:
    solver<SType,Type,MeshType>::directSolver(dict,fvMsh,l),
    APtrs_(MeshType::numberOfDirections),
    pivotPtrs_(MeshType::numberOfDirections),
    gCellNumbers_(fvMsh,l)
{}

template<class SType, class Type, class MeshType>
void APLU<SType,Type,MeshType>::prepare
(
    linearSystem<SType,Type,MeshType>& xEqn
)
{
    const globalRowData<SType,Type,MeshType> gRowData(xEqn, this->l_);

    const List<bool> singular(xEqn.singular());
    const List<bool> diagonal(xEqn.diagonal());

    if (Pstream::master())
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            // Nothing to be done for diagonal systems

            if (diagonal[d])
                continue;

            const auto& cellNumbers = gCellNumbers_.data(d);
            const auto& rowData = gRowData.data(d);

            const label n = gCellNumbers_.numberOfCells(d);

            const bool symm =
                SType::csType::typeName == symmStencil::csType::typeName;

            // Create linear system

            APtrs_.set(d, new scalarSquareMatrix(n));
            scalarSquareMatrix& A = APtrs_[d];
            A = Zero;

            label offset = 0;
            forAll(cellNumbers, proc)
            {
                const auto& nums = cellNumbers[proc];
                const auto& rows = rowData[proc];

                forAll(nums, cell)
                {
                    const label i = offset+cell;

                    const SType S = rows[cell];

                    for (int s = 0; s < SType::nComponents; s++)
                    {
                        const label j = nums[cell][s];

                        if (j > -1)
                        {
                            A[i][j] = S[s];

                            if (symm && i != j)
                                A[j][i] = S[s];
                        }
                        else if (Foam::mag(S[s]) > 1e-8)
                        {
                            // Homogeneous Neumann in case of non-eliminated
                            // boundary coefficient

                            A[i][i] += S[s];
                        }
                    }
                }

                offset += nums.size();
            }

            // For a singular matrix force the first value to zero

            if (singular[d])
            {
                for (int i = 0; i < n; i++)
                    A[i][0] = A[0][i] = 0.0;

                A[0][0] = 1.0;
            }

            // Print the full matrix for debugging

            // for (int i = 0; i < n; i++)
            // {
            //     for (int j = 0; j < n; j++)
            //         Info<< A(i,j) << " ";
            //     Info<< endl;
            // }

            // Decompose matrix

            pivotPtrs_.set(d, new labelList(n));
            labelList& pivots = pivotPtrs_[d];

            LUDecompose(A, pivots);
        }
    }

    solver<SType,Type,MeshType>::directSolver::prepare(xEqn);
}

template<class SType, class Type, class MeshType>
void APLU<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& xEqn
)
{
    const globalRowData<SType,Type,MeshType> gRowData(xEqn, this->l_);

    meshLevel<Type,MeshType>& x = xEqn.x()[this->l_];
    const meshLevel<Type,MeshType>& b = xEqn.b()[this->l_];
    const meshLevel<SType,MeshType>& A = xEqn.A()[this->l_];

    const List<bool> singular(xEqn.singular());
    const List<bool> diagonal(xEqn.diagonal());

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        if (diagonal[d])
        {
            forAllCells(x[d], i, j, k)
                x[d](i,j,k) = b[d](i,j,k)/A[d](i,j,k).center();
        }
        else if (Pstream::master())
        {
            const auto& cellNumbers = gCellNumbers_.data(d);
            const label n = gCellNumbers_.numberOfCells(d);

            // Get the right-hand side

            List<Type> rhs(n);

            label offset = 0;
            forAll(cellNumbers, proc)
            {
                const auto& nums = cellNumbers[proc];

                if (proc == Pstream::myProcNo())
                {
                    int c = 0;
                    forAllCells(x[d], i, j, k)
                        rhs[offset + c++] = b[d](i,j,k);
                }
                else
                {
                    UIPstream::read
                    (
                        Pstream::commsTypes::blocking,
                        proc,
                        reinterpret_cast<char*>(&rhs[offset]),
                        nums.size()*sizeof(Type),
                        0,
                        UPstream::worldComm
                    );
                }

                offset += nums.size();
            }

            // For a singular matrix the first value is forced to zero

            if (singular[d])
                rhs[0] = Zero;

            // Solve

            const scalarSquareMatrix& A = APtrs_[d];
            const labelList& pivots = pivotPtrs_[d];

            LUBacksubstitute(A, pivots, rhs);

            // Set mean to zero

            if (singular[d])
            {
                Type avg = average(rhs);

                forAll(rhs, i)
                    rhs[i] -= avg;
            }

            // Send/copy solutions

            offset = 0;
            forAll(cellNumbers, proc)
            {
                const auto& nums = cellNumbers[proc];

                if (proc != Pstream::myProcNo())
                {
                    UOPstream::write
                    (
                        Pstream::commsTypes::blocking,
                        proc,
                        reinterpret_cast<char*>(&rhs[offset]),
                        nums.size()*sizeof(Type),
                        0,
                        UPstream::worldComm
                    );
                }
                else
                {
                    label c = 0;
                    forAllCells(x[d], i, j, k)
                        x[d](i,j,k) = rhs[offset + c++];
                }

                offset += nums.size();
            }
        }
        else
        {
            // Send right-hand side

            List<Type> rhs(cmptProduct(x[d].N()));

            int c = 0;
            forAllCells(x[d], i, j, k)
                rhs[c++] = b[d](i,j,k);

            UOPstream::write
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo(),
                reinterpret_cast<char*>(rhs.begin()),
                rhs.byteSize(),
                0,
                UPstream::worldComm
            );

            // Receive and copy solution

            UIPstream::read
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo(),
                reinterpret_cast<char*>(rhs.begin()),
                rhs.byteSize(),
                0,
                UPstream::worldComm
            );

            c = 0;
            forAllCells(x[d], i, j, k)
                x[d](i,j,k) = rhs[c++];
        }
    }

    x.correctBoundaryConditions();
}

}

}

}
