#include "Eigen.H"

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
Eigen<SType,Type,MeshType>::Eigen
(
    const dictionary& dict,
    const fvMesh& fvMsh,
    const label l
)
:
    solver<SType,Type,MeshType>::directSolver(dict,fvMsh,l),
    gCellNumbers_(fvMsh,l),
    APtrs_(MeshType::numberOfDirections),
    solverPtrs_(MeshType::numberOfDirections)
{
    if (Pstream::master())
        for (int d = 0; d < MeshType::numberOfDirections; d++)
            solverPtrs_.set
            (
                d,
                EigenSolver::New
                (
                    dict.lookupOrDefault<word>("EigenSolver", "default")
                ).ptr()
            );
}

template<class SType, class Type, class MeshType>
void Eigen<SType,Type,MeshType>::prepare
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

            bool symm =
                SType::csType::typeName == symmStencil::csType::typeName;

            // Create linear system. Increase dimension by one, to allow for
            // singular system augmentation.

            typedef ::Eigen::Triplet<double> tType;

            std::vector<tType> coeffs;
            coeffs.reserve(n*(SType::nComponents + singular[d]*2) + 1);

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
                            coeffs.push_back(tType(i,j,S[s]));

                            if (symm && i != j)
                                coeffs.push_back(tType(j,i,S[s]));
                        }
                        else if (Foam::mag(S[s]) > 1e-8)
                        {
                            // Homogeneous Neumann in case of non-eliminated
                            // boundary coefficient

                            coeffs.push_back(tType(i,i,S[s]));
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
                    coeffs.push_back(tType(i,n,1));
                    coeffs.push_back(tType(n,i,1));
                }
            }
            else
            {
                // Set the auxilary variable to zero

                coeffs.push_back(tType(n,n,1));
            }

            // Set matrix and compute decomposition

            APtrs_.set(d, new EigenSolver::matrixType(n+1,n+1));
            EigenSolver::matrixType& A = APtrs_[d];

            A.setFromTriplets(coeffs.begin(), coeffs.end());
            A.makeCompressed();

            solverPtrs_[d].prepare(A);
        }
    }

    solver<SType,Type,MeshType>::directSolver::prepare(xEqn, singular);
}

template<class SType, class Type, class MeshType>
void Eigen<SType,Type,MeshType>::solve
(
    linearSystem<SType,Type,MeshType>& xEqn,
    const List<bool>& singular
)
{
    const globalRowData<SType,Type,MeshType> gRowData(xEqn, this->l_);

    meshLevel<Type,MeshType>& x = xEqn.x()[this->l_];
    const meshLevel<Type,MeshType>& b = xEqn.b()[this->l_];

    const label m = list(pTraits<Type>::one).size();

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

            // Copy to Eigen right-hand side type

            EigenSolver::rhsType B(n+1,m);

            forAll(rhs, i)
            {
                const scalarList L(list(rhs[i]));

                forAll(L, j)
                    B(i,j) = L[j];
            }

            // Solve

            EigenSolver::rhsType X;
            solverPtrs_[d].solve(X,B);

            // Copy back to rhs

            forAll(rhs, i)
            {
                scalarList L(m);

                forAll(L, j)
                    L[j] = X(i,j);

                rhs[i] = listToType<Type>(L);
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
