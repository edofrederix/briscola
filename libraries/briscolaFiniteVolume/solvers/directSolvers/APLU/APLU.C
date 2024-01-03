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
    const linearSystem<SType,Type,MeshType>& xEqn,
    const List<bool>& singular
)
{
    const globalRowData<SType,Type,MeshType> gRowData(xEqn, this->l_);

    if (Pstream::master())
    {
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const auto& cellNumbers = gCellNumbers_.data(d);
            const auto& rowData = gRowData.data(d);

            const label n = gCellNumbers_.numberOfCells(d);

            const bool symm =
                SType::csType::typeName == symmStencil::csType::typeName;

            // Create linear system. Increase dimension by one, to allow for
            // singular system augmentation.

            APtrs_.set(d, new scalarSquareMatrix(n+1));
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

            if (singular[d])
            {
                // Singular matrix, augment the system (see "Multigrid", U.
                // Trottenberg et al., 2001)

                for (int i = 0; i < n; i++)
                {
                    A(i,n) = 1.0;
                    A(n,i) = 1.0;
                }
            }
            else
            {
                // Set the auxilary variable to zero

                A(n,n) = 1.0;
            }

            // Decompose matrix

            pivotPtrs_.set(d, new labelList(n+1));
            labelList& pivots = pivotPtrs_[d];

            LUDecompose(A, pivots);
        }
    }

    solver<SType,Type,MeshType>::directSolver::prepare(xEqn, singular);
}

template<class SType, class Type, class MeshType>
void APLU<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& xEqn,
    const List<bool>& singular
)
{
    const globalRowData<SType,Type,MeshType> gRowData(xEqn, this->l_);

    meshLevel<Type,MeshType>& x = xEqn.x()[this->l_];
    const meshLevel<Type,MeshType>& b = xEqn.b()[this->l_];

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        if (Pstream::master())
        {
            const auto& cellNumbers = gCellNumbers_.data(d);
            const label n = gCellNumbers_.numberOfCells(d);

            // Get the right-hand side

            List<Type> rhs(n+1);

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

            rhs[n] = Zero;

            // Solve

            const scalarSquareMatrix& A = APtrs_[d];
            const labelList& pivots = pivotPtrs_[d];

            LUBacksubstitute(A, pivots, rhs);

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
