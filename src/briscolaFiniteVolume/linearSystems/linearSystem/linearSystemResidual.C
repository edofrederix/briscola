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
void linearSystem<SType,Type,MeshType>::Residual
(
    meshDirection<Type,MeshType>& res,
    const faceLabel& paddingMask,
    const label padding
) const
{
    setForcingMask();

    const label l = res.levelNum();
    const label d = res.directionNum();

    block<Type>& B = res.B();
    const labelVector shape = B.shape();

    // Restricted array pointers

    Type* const __restrict__ res_arr = B.begin();

    const Type* const __restrict__ x_arr = this->x()[l][d].B().begin();
    const Type* const __restrict__ b_arr = this->b()[l][d].B().begin();
    const label* const __restrict__ f_arr =
        this->forcingMask()[l][d].B().begin();

    // Reinterpret the matrix as a contiguous array of scalars

    const scalar* const __restrict__ A_arr =
        reinterpret_cast<const scalar*>(this->A()[l][d].B().begin());

    const faceLabel I = res.I();

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

    // Compute residuals

    int c = lin(S, shape);

    if (!Mask || !this->x().immersedBoundaryConditions().size())
    {
        // Compute residual for each cell

        for (int i = 0; i < M.x(); i++)
        {
            for (int j = 0; j < M.y(); j++)
            {
                for (int k = 0; k < M.z(); k++)
                {
                    res_arr[c] =
                        b_arr[c]
                      - P(A_arr, x_arr, c, S_i, S_j, S_k);

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
        // Compute residual for each unmasked cell

        for (int i = 0; i < M.x(); i++)
        {
            for (int j = 0; j < M.y(); j++)
            {
                for (int k = 0; k < M.z(); k++)
                {
                    if (f_arr[c])
                    {
                        res_arr[c] = Zero;
                    }
                    else
                    {
                        res_arr[c] =
                            b_arr[c]
                          - P(A_arr, x_arr, c, S_i, S_j, S_k);
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
void linearSystem<SType,Type,MeshType>::ResidualBoundary
(
    meshDirection<Type,MeshType>& res,
    const faceLabel& paddingMask,
    const label padding
) const
{
    setForcingMask();

    const label l = res.levelNum();
    const label d = res.directionNum();

    block<Type>& B = res.B();
    const labelVector shape = B.shape();

    // Restricted array pointers

    Type* const __restrict__ res_arr = B.begin();

    const Type* const __restrict__ x_arr = this->x()[l][d].B().begin();
    const Type* const __restrict__ b_arr = this->b()[l][d].B().begin();
    const label* const __restrict__ f_arr =
        this->forcingMask()[l][d].B().begin();

    // Reinterpret the matrix as a contiguous array of scalars

    const scalar* const __restrict__ A_arr =
        reinterpret_cast<const scalar*>(this->A()[l][d].B().begin());

    const faceLabel I = res.I();

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

        // Compute residuals

        int c = lin(S, shape);

        if (!Mask || !this->x().immersedBoundaryConditions().size())
        {
            // Compute residual for each cell

            for (int i = 0; i < M.x(); i++)
            {
                for (int j = 0; j < M.y(); j++)
                {
                    for (int k = 0; k < M.z(); k++)
                    {
                        res_arr[c] =
                            b_arr[c]
                          - P(A_arr, x_arr, c, S_i, S_j, S_k);

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
            // Compute residual for each unmasked cell

            for (int i = 0; i < M.x(); i++)
            {
                for (int j = 0; j < M.y(); j++)
                {
                    for (int k = 0; k < M.z(); k++)
                    {
                        if (f_arr[c])
                        {
                            res_arr[c] = Zero;
                        }
                        else
                        {
                            res_arr[c] =
                                b_arr[c]
                              - P(A_arr, x_arr, c, S_i, S_j, S_k);
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
void linearSystem<SType,Type,MeshType>::residual
(
    meshField<Type, MeshType>& res
) const
{
    forAll(res, l)
        this->residual<CorrBcs,Mask>(res[l]);
}

template<class SType, class Type, class MeshType>
template<bool CorrBcs, bool Mask>
void linearSystem<SType,Type,MeshType>::residual
(
    meshLevel<Type,MeshType>& res
) const
{
    meshLevel<Type,MeshType>& x =
        const_cast<meshLevel<Type,MeshType>&>(this->x()[res.levelNum()]);

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

        // Compute residual for internal cells

        forAll(x, d)
            Residual(res[d], paddingMask, 1);

        // Evaluate boundary conditions

        if (!eliminated_)
            x.template evaluate<eliminatedBcs,1>();

        if (Pstream::parRun())
            UPstream::waitRequests(nReq);

        x.template evaluate<nonEliminatedBcs,1>();

        // Compute residual for boundary cells

        forAll(x, d)
            ResidualBoundary(res[d], paddingMask, 1);
    }
    else
    {
        forAll(x, d)
            Residual(res[d], paddingMask);
    }
}

template<class SType, class Type, class MeshType>
template<bool CorrBcs, bool Mask>
tmp<meshField<Type, MeshType>>
linearSystem<SType,Type,MeshType>::residual() const
{
    tmp<meshField<Type,MeshType>> tRes =
        meshField<Type,MeshType>::New
        (
            "residual",
            fvMsh_
        );

    this->residual<CorrBcs,Mask>(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
template<bool CorrBcs, bool Mask>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::residual(const label l) const
{
    tmp<meshLevel<Type,MeshType>> tRes =
        meshLevel<Type,MeshType>::New(fvMsh_,l);

    this->residual<CorrBcs,Mask>(tRes.ref());

    return tRes;
}

}

}

}
