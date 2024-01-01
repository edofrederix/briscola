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
    gCellNumbers_(fvMsh,l)
{}

template<class SType, class Type, class MeshType>
void APLU<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& xEqn,
    const List<bool>& singular
)
{
    const globalRowData<SType,Type,MeshType> gRowData(xEqn, this->l_);

    meshLevel<Type,MeshType>& x = xEqn.x()[this->l_];

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        if (Pstream::master())
        {
            const auto& cellNumbers = gCellNumbers_.data(d);
            const auto& rowData = gRowData.data(d);

            const label n = gCellNumbers_.numberOfCells(d);

            // Create linear system. Increase dimension by one, to allow for
            // singular system augmentation.

            scalarSquareMatrix A(n+1, Zero);
            List<Type> b(n+1);
            b[n] = Zero;

            const bool symm =
                SType::csType::typeName == symmStencil::csType::typeName;

            label offset = 0;
            forAll(cellNumbers, proc)
            {
                const auto& nums = cellNumbers[proc];
                const auto& rows = rowData[proc];

                forAll(nums, cell)
                {
                    b[offset+cell] = rows[cell].source();

                    const SType S = rows[cell].stencil();

                    for (int s = 0; s < SType::nComponents; s++)
                    {
                        if (nums[cell][s] > -1)
                        {
                            A[offset+cell][nums[cell][s]] = S[s];

                            if (symm)
                                A[nums[cell][s]][offset+cell] =
                                    A[offset+cell][nums[cell][s]];
                        }
                        else if (Foam::mag(S[s]) > 1e-8)
                        {
                            // Homogeneous Neumann in case of non-eliminated
                            // boundary coefficient

                            A[offset+cell][nums[cell][s]] += S[s];
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

            // Solve

            LUsolve(A,b);

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
                        reinterpret_cast<char*>(&b[offset]),
                        nums.size()*sizeof(Type),
                        0,
                        UPstream::worldComm
                    );
                }
                else
                {
                    label c = 0;
                    forAllCells(x[d], i, j, k)
                        x[d](i,j,k) = b[offset + c++];
                }

                offset += nums.size();
            }
        }
        else
        {
            // Receive and copy solution

            List<Type> b(cmptProduct(x[d].N()));

            UIPstream::read
            (
                Pstream::commsTypes::blocking,
                Pstream::masterNo(),
                reinterpret_cast<char*>(b.begin()),
                b.byteSize(),
                0,
                UPstream::worldComm
            );

            label c = 0;
            forAllCells(x[d], i, j, k)
                x[d](i,j,k) = b[c++];
        }
    }

    x.correctBoundaryConditions();
}

}

}

}
