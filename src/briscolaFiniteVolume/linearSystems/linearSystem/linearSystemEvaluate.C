#include "linearSystem.H"
#include "linearSystemFunctions.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
template<bool Mask>
void linearSystem<SType,Type,MeshType>::Evaluate
(
    meshDirection<Type,MeshType>& eval,
    const faceLabel& paddingMask,
    const label padding
) const
{
    setForcingMask();

    const label l = eval.levelNum();
    const label d = eval.directionNum();

    block<Type>& B = eval.B();
    const labelVector shape = B.shape();

    // Restricted array pointers

    Type* const __restrict__ eval_arr = B.begin();

    const Type* const __restrict__ x_arr = this->x()[l][d].B().begin();
    const Type* const __restrict__ b_arr = this->b()[l][d].B().begin();
    const label* const __restrict__ f_arr =
        this->forcingMask()[l][d].B().begin();

    const scalar* const __restrict__ icv_arr =
        this->fvMsh_.template
        metrics<MeshType>().inverseCellVolumes()[l][d].B().begin();

    // Reinterpret the matrix as a contiguous array of scalars

    const scalar* const __restrict__ A_arr =
        reinterpret_cast<const scalar*>(this->A()[l][d].B().begin());

    const faceLabel I = eval.I();

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

    // Row product function

    linearSystemFun::rowProduct<SType,Type> P;

    // Compute evaluations

    int c = lin(S, shape);

    if (!Mask || !this->x().immersedBoundaryConditions().size())
    {
        // Compute evaluation for each cell

        for (int i = 0; i < M.x(); i++)
        {
            for (int j = 0; j < M.y(); j++)
            {
                for (int k = 0; k < M.z(); k++)
                {
                    eval_arr[c] =
                        P(A_arr, x_arr, c, S_i, S_j, S_k)
                      - b_arr[c];

                    eval_arr[c] *= icv_arr[c];

                    // Jump to next cell
                    c++;
                }

                // Jump to next line
                c += J_k;
            }

            // Jump to next plane
            c += J_j;
        }
    }
    else
    {
        // Compute evaluation for each unmasked cell

        for (int i = 0; i < M.x(); i++)
        {
            for (int j = 0; j < M.y(); j++)
            {
                for (int k = 0; k < M.z(); k++)
                {
                    if (f_arr[c])
                    {
                        eval_arr[c] = Zero;
                    }
                    else
                    {
                        eval_arr[c] =
                            P(A_arr, x_arr, c, S_i, S_j, S_k)
                          - b_arr[c];

                        eval_arr[c] *= icv_arr[c];
                    }

                    // Jump to next cell
                    c++;
                }

                // Jump to next line
                c += J_k;
            }

            // Jump to next plane
            c += J_j;
        }
    }
}

template<class SType, class Type, class MeshType>
template<bool Mask>
void linearSystem<SType,Type,MeshType>::EvaluateBoundary
(
    meshDirection<Type,MeshType>& eval,
    const faceLabel& paddingMask,
    const label padding
) const
{
    setForcingMask();

    const label l = eval.levelNum();
    const label d = eval.directionNum();

    block<Type>& B = eval.B();
    const labelVector shape = B.shape();

    // Restricted array pointers

    Type* const __restrict__ eval_arr = B.begin();

    const Type* const __restrict__ x_arr = this->x()[l][d].B().begin();
    const Type* const __restrict__ b_arr = this->b()[l][d].B().begin();
    const label* const __restrict__ f_arr =
        this->forcingMask()[l][d].B().begin();

    const scalar* const __restrict__ icv_arr =
        this->fvMsh_.template
        metrics<MeshType>().inverseCellVolumes()[l][d].B().begin();

    // Reinterpret the matrix as a contiguous array of scalars

    const scalar* const __restrict__ A_arr =
        reinterpret_cast<const scalar*>(this->A()[l][d].B().begin());

    const faceLabel I = eval.I();

    // Row product function

    linearSystemFun::rowProduct<SType,Type> P;

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

        // Compute evaluations

        int c = lin(S, shape);

        if (!Mask || !this->x().immersedBoundaryConditions().size())
        {
            // Compute evaluation for each cell

            for (int i = 0; i < M.x(); i++)
            {
                for (int j = 0; j < M.y(); j++)
                {
                    for (int k = 0; k < M.z(); k++)
                    {
                        eval_arr[c] =
                            P(A_arr, x_arr, c, S_i, S_j, S_k)
                          - b_arr[c];

                        eval_arr[c] *= icv_arr[c];

                        // Jump to next cell
                        c++;
                    }

                    // Jump to next line
                    c += J_k;
                }

                // Jump to next plane
                c += J_j;
            }
        }
        else
        {
            // Compute evaluation for each unmasked cell

            for (int i = 0; i < M.x(); i++)
            {
                for (int j = 0; j < M.y(); j++)
                {
                    for (int k = 0; k < M.z(); k++)
                    {
                        if (f_arr[c])
                        {
                            eval_arr[c] = Zero;
                        }
                        else
                        {
                            eval_arr[c] =
                                P(A_arr, x_arr, c, S_i, S_j, S_k)
                              - b_arr[c];

                            eval_arr[c] *= icv_arr[c];
                        }

                        // Jump to next cell
                        c++;
                    }

                    // Jump to next line
                    c += J_k;
                }

                // Jump to next plane
                c += J_j;
            }
        }
    }
}

template<class SType, class Type, class MeshType>
template<bool CorrBcs, bool Mask>
void linearSystem<SType,Type,MeshType>::evaluate
(
    meshField<Type, MeshType>& eval
) const
{
    forAll(eval, l)
        this->evaluate<CorrBcs,Mask>(eval[l]);
}

template<class SType, class Type, class MeshType>
template<bool CorrBcs, bool Mask>
void linearSystem<SType,Type,MeshType>::evaluate
(
    meshLevel<Type,MeshType>& eval
) const
{
    meshLevel<Type,MeshType>& x =
        const_cast<meshLevel<Type,MeshType>&>(this->x()[eval.levelNum()]);

    linearSystem<SType,Type,MeshType>& sys =
        const_cast<linearSystem<SType,Type,MeshType>&>(*this);

    const faceLabel paddingMask =
        faceLabel::one - sys.eliminatedBoundaryMasks()[x.levelNum()];

    if (CorrBcs)
    {
        // Prepare boundary condition transfers

        const label nReq = Pstream::nRequests();

        x.template prepare<nonEliminatedBcs,1>();

        if (!eliminated_)
            x.template prepare<eliminatedBcs,1>();

        // Evaluate internal cells

        forAll(x, d)
            Evaluate(eval[d], paddingMask, 1);

        // Evaluate boundary conditions

        if (!eliminated_)
            x.template evaluate<eliminatedBcs,1>();

        if (Pstream::parRun())
            UPstream::waitRequests(nReq);

        x.template evaluate<nonEliminatedBcs,1>();

        // Evaluate boundary cells

        forAll(x, d)
            EvaluateBoundary(eval[d], paddingMask, 1);
    }
    else
    {
        forAll(x, d)
            Evaluate(eval[d], paddingMask);
    }
}

template<class SType, class Type, class MeshType>
template<bool CorrBcs, bool Mask>
tmp<meshField<Type, MeshType>>
linearSystem<SType,Type,MeshType>::evaluate() const
{
    tmp<meshField<Type,MeshType>> tEval =
        meshField<Type,MeshType>::New
        (
            "evaluate",
            fvMsh_
        );

    this->evaluate<CorrBcs,Mask>(tEval.ref());

    return tEval;
}

template<class SType, class Type, class MeshType>
template<bool CorrBcs, bool Mask>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::evaluate(const label l) const
{
    tmp<meshLevel<Type,MeshType>> tEval =
        meshLevel<Type,MeshType>::New(fvMsh_,l);

    this->evaluate<CorrBcs,Mask>(tEval.ref());

    return tEval;
}

}

}

}
