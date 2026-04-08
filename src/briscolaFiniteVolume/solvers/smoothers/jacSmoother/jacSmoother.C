#include "jacSmoother.H"
#include "linearSystemFunctions.H"
#include "diagonalSmoother.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
const scalar jacSmoother<SType,Type,MeshType>::jacSmoother::omega_ = 0.8;

template<class SType, class Type, class MeshType>
inline void jacSmoother<SType,Type,MeshType>::jacSmoother::Sweep
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label d
)
{
    if (sys.fvMsh()[l].empty())
        return;

    block<Type>& B = sys.x()[l][d].B();
    const labelVector shape = B.shape();

    // Restricted array pointers

    Type* const __restrict__ x_arr = B.begin();

    const Type* const __restrict__ b_arr = sys.b()[l][d].B().begin();
    const label* const __restrict__ f_arr = sys.forcingMask()[l][d].B().begin();

    // Reinterpret the matrix as a contiguous array of scalars

    const scalar* const __restrict__ A_arr =
        reinterpret_cast<const scalar*>(sys.A()[l][d].B().begin());

    const faceLabel I = sys.x()[l][d].I();

    // Data in B is padded by ghosts

    const labelVector S = I.lower() + G*unitXYZ;
    const labelVector E = I.upper() + G*unitXYZ;
    const labelVector M = E - S;

    // Array to temporarily store solution

    block<Type> y(M);

    Type* const __restrict__ y_arr = y.begin();

    // Strides in i and j (data is contiguous in k)

    const label S_i = lin(S+unitX, shape) - lin(S, shape);
    const label S_j = lin(S+unitY, shape) - lin(S, shape);
    const label S_k = 1;

    // Jump after each line in k and plane in (j,k)

    const label J_k = lin(S+unitY, shape) - lin(S+unitZ*M.z(), shape);
    const label J_j = lin(S+unitX, shape) - lin(S+unitY*M.y(), shape);

    // Off-diagonal row product function

    linearSystemFun::offDiagRowProduct<SType,Type> P;

    // Sweep

    int c = lin(S, shape);
    int e = 0;

    for (int i = 0; i < M.x(); i++)
    {
        for (int j = 0; j < M.y(); j++)
        {
            for (int k = 0; k < M.z(); k++)
            {
                if (!f_arr[c])
                    y_arr[e] =
                        (1.0 - omega_)*x_arr[c]
                      + omega_
                      * (
                            b_arr[c]
                          - P(A_arr, x_arr, c, S_i, S_j, S_k)
                        )
                      / A_arr[c*N];

                // Jump to next cell
                c++;
                e++;
            }

            // Jump to next line
            c += J_k;
        }

        // Jump to next plane
        c += J_j;
    }

    // Copy back

    c = lin(S, shape);
    e = 0;

    for (int i = 0; i < M.x(); i++)
    {
        for (int j = 0; j < M.y(); j++)
        {
            for (int k = 0; k < M.z(); k++)
                x_arr[c++] = y_arr[e++];

            // Jump to next line
            c += J_k;
        }

        // Jump to next plane
        c += J_j;
    }
}

template<class SType, class Type, class MeshType>
void jacSmoother<SType,Type,MeshType>::jacSmoother::Smooth
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label sweeps,
    const labelList& converged
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        if (l==0)
            x.correctImmersedBoundaryConditions();

        x.template correct<nonEliminatedBcs>();

        if (!sys.eliminated())
            x.template correct<eliminatedBcs>();

        forAll(x, d)
        if (!converged[d])
        {
            if (sys.diagonal())
            {
                diagonalSmoother<SType,Type,MeshType>::Sweep(sys, l, d);
            }
            else
            {
                Sweep(sys, l, d);
            }
        }
    }
}

template<class SType, class Type, class MeshType>
jacSmoother<SType,Type,MeshType>::jacSmoother
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>::smoother(dict,fvMsh)
{}

template<class SType, class Type, class MeshType>
void jacSmoother<SType,Type,MeshType>::jacSmoother::smooth
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
