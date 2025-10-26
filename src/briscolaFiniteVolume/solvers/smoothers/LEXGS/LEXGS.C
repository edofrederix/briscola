#include "LEXGS.H"
#include "linearSystemFunctions.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
inline void LEXGS<SType,Type,MeshType>::LEXGS::Sweep
(
    meshDirection<Type,MeshType>& x,
    const meshDirection<SType,MeshType>& A,
    const meshDirection<Type,MeshType>& b,
    const meshDirection<label,MeshType>& f
)
{
    block<Type>& B = x.B();
    const labelVector shape = B.shape();

    // Restricted array pointers

    Type* const __restrict__ x_arr = B.begin();

    const Type* const __restrict__ b_arr = b.B().begin();
    const label* const __restrict__ f_arr = f.B().begin();

    // Reinterpret the matrix as a contiguous array of scalars

    const scalar* const __restrict__ A_arr =
        reinterpret_cast<const scalar*>(A.B().begin());

    const faceLabel I = x.I();

    // Data in B is padded by ghosts

    const labelVector S = I.lower() + G*unitXYZ;
    const labelVector E = I.upper() + G*unitXYZ;
    const labelVector M = E - S;

    // Strides in i and j (data is contiguous in k)

    const label S_i = lin(S+unitX, shape) - lin(S, shape);
    const label S_j = lin(S+unitY, shape) - lin(S, shape);
    const label S_k = 1;

    // Jump after each line in k and plane in (j,k)

    const label J_k = lin(S+unitY, shape) - lin(S+unitZ*M.z(), shape);
    const label J_j = lin(S+unitX, shape) - lin(S+unitY*M.y(), shape);

    // Off-diagonal row product function

    linearSystemFun::offDiagRowProduct<SType,Type> P;

    // Forward sweep

    int c = lin(S, shape);

    for (int i = 0; i < M.x(); i++)
    {
        for (int j = 0; j < M.y(); j++)
        {
            for (int k = 0; k < M.z(); k++)
            {
                if (!f_arr[c])
                    x_arr[c] =
                        (
                            b_arr[c]
                          - P(A_arr, x_arr, c, S_i, S_j, S_k)
                        )
                      / A_arr[c*N];

                // Jump to next cell
                c++;
            }

            // Jump to next line
            c += J_k;
        }

        // Jump to next plane
        c += J_j;
    }

    // Backward sweep

    c = lin(E-unitXYZ, shape);

    for (int i = 0; i < M.x(); i++)
    {
        for (int j = 0; j < M.y(); j++)
        {
            for (int k = 0; k < M.z(); k++)
            {
                if (!f_arr[c])
                    x_arr[c] =
                        (
                            b_arr[c]
                          - P(A_arr, x_arr, c, S_i, S_j, S_k)
                        )
                      / A_arr[c*N];

                // Jump to next cell
                c--;
            }

            // Jump to next line
            c -= J_k;
        }

        // Jump to next plane
        c -= J_j;
    }
}

template<class SType, class Type, class MeshType>
void LEXGS<SType,Type,MeshType>::LEXGS::Smooth
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label sweeps,
    const labelList& converged
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const meshLevel<SType,MeshType>& A = sys.A()[l];
    const meshLevel<Type,MeshType>& b = sys.b()[l];
    const meshLevel<label,MeshType>& f = sys.forcingMask()[l];

    const List<bool> diagonal(sys.diagonal());

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        forAll(x, d)
        if (!converged[d])
        {
            if (diagonal[d])
            {
                solver<SType,Type,MeshType>::smoother::smoothDiag(sys, l, d);
            }
            else
            {
                Sweep(x[d], A[d], b[d], f[d]);
            }
        }

        if (l == 0)
            x.correctImmersedBoundaryConditions();

        x.correctNonEliminatedBoundaryConditions();
    }

    x.correctEliminatedBoundaryConditions();
    x.correctUnsetBoundaryConditions();
}

template<class SType, class Type, class MeshType>
LEXGS<SType,Type,MeshType>::LEXGS
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>::smoother(dict,fvMsh)
{}

template<class SType, class Type, class MeshType>
void LEXGS<SType,Type,MeshType>::LEXGS::smooth
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label sweeps,
    const labelList& converged
)
{
    this->Smooth(sys, l, sweeps, converged);
}

}

}

}
