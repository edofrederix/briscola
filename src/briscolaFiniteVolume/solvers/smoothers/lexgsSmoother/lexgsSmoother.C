#include "lexgsSmoother.H"
#include "linearSystemFunctions.H"
#include "diagonalSmoother.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
inline void lexgsSmoother<SType,Type,MeshType>::lexgsSmoother::Sweep
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label d,
    const faceLabel& paddingMask,
    const label padding
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

    // Data in B is padded by ghosts and by the boundary mask

    const labelVector S = I.lower() + G*unitXYZ + padding*paddingMask.lower();
    const labelVector E = I.upper() + G*unitXYZ - padding*paddingMask.upper();
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
inline void lexgsSmoother<SType,Type,MeshType>::lexgsSmoother::SweepBoundary
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label d,
    const faceLabel& paddingMask,
    const label padding
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

    // Off-diagonal row product function

    linearSystemFun::offDiagRowProduct<SType,Type> P;

    for (int f = 0; f < 6; f++)
    if (paddingMask[f])
    {
        const label fd = f/2;

        // Data in B is padded by ghosts

        labelVector S = I.lower() + G*unitXYZ;
        labelVector E = I.upper() + G*unitXYZ;

        // Adjust upper or lower limits

        if (f%2 == 0)
        {
            E[fd] = S[fd] + padding;
        }
        else
        {
            S[fd] = E[fd] - padding;
        }

        const labelVector M = E - S;

        // Strides in i and j (data is contiguous in k)

        const label S_i = lin(S+unitX, shape) - lin(S, shape);
        const label S_j = lin(S+unitY, shape) - lin(S, shape);
        const label S_k = 1;

        // Jump after each line in k and plane in (j,k)

        const label J_k = lin(S+unitY, shape) - lin(S+unitZ*M.z(), shape);
        const label J_j = lin(S+unitX, shape) - lin(S+unitY*M.y(), shape);

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
}

template<class SType, class Type, class MeshType>
void lexgsSmoother<SType,Type,MeshType>::lexgsSmoother::Smooth
(
    linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label sweeps,
    const labelList& converged
)
{
    meshLevel<Type,MeshType>& x = sys.x()[l];

    const faceLabel paddingMask =
        faceLabel::one - sys.eliminatedBoundaryMasks()[l];

    for (label sweep = 0; sweep < sweeps; sweep++)
    {
        if (l == 0)
            x.correctImmersedBoundaryConditions();

        if (sys.diagonal())
        {
            // No boundary condition correction needed

            forAll(x, d)
                if (!converged[d])
                    diagonalSmoother<SType,Type,MeshType>::Sweep(sys, l, d);
        }
        else
        {
            // Prepare boundary condition transfers

            const label nReq = Pstream::nRequests();

            x.template prepare<nonEliminatedBcs,1>();

            if (!sys.eliminated())
                x.template prepare<eliminatedBcs,1>();

            // Sweep internal cells

            forAll(x, d)
                if (!converged[d])
                    Sweep(sys, l, d, paddingMask, 1);

            // Evaluate boundary conditions

            if (!sys.eliminated())
                x.template evaluate<eliminatedBcs,1>();

            if (Pstream::parRun())
                UPstream::waitRequests(nReq);

            x.template evaluate<nonEliminatedBcs,1>();

            // Sweep boundary cells

            forAll(x, d)
                if (!converged[d])
                    SweepBoundary(sys, l, d, paddingMask, 1);
        }
    }
}

template<class SType, class Type, class MeshType>
lexgsSmoother<SType,Type,MeshType>::lexgsSmoother
(
    const dictionary& dict,
    const fvMesh& fvMsh
)
:
    solver<SType,Type,MeshType>::smoother(dict,fvMsh)
{}

template<class SType, class Type, class MeshType>
void lexgsSmoother<SType,Type,MeshType>::lexgsSmoother::smooth
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
