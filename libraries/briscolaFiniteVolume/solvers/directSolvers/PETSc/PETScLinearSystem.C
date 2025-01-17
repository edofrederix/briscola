#include "PETScLinearSystem.H"
#include "PstreamGlobalsLsa.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
PETScLinearSystem<SType,Type,MeshType>::PETScLinearSystem
(
    const linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label nParts
)
:
    linearSystemAggregation<SType,Type,MeshType>(sys, l, nParts),
    sys_(sys),
    matrices_(MeshType::numberOfDirections),
    l_(l)
{}

template<class SType, class Type, class MeshType>
PETScLinearSystem<SType,Type,MeshType>::PETScLinearSystem
(
    const PETScLinearSystem<SType,Type,MeshType>& A
)
:
    linearSystemAggregation<SType,Type,MeshType>(A),
    sys_(A.sys_),
    matrices_(A.matrices_),
    l_(A.l_)
{}

template<class SType, class Type, class MeshType>
PETScLinearSystem<SType,Type,MeshType>::~PETScLinearSystem()
{}

template<class SType, class Type, class MeshType>
void PETScLinearSystem<SType,Type,MeshType>::prepare(const label d)
{
    typedef typename SType::fullStencilType FullSType;

    const linearSystemAggregation<SType,Type,MeshType>& lsa = *this;

    const MPI_Comm& comm =
        Pstream::parRun()
      ? PstreamGlobals::lsaGetComm(lsa.masterCommNum())
      : PETSC_COMM_SELF;

    List<List<FullSType>> coeffs;
    lsa.rowCoeffs(coeffs, sys_, d);

    List<List<FixedList<label,FullSType::nComponents>>> colNums =
        lsa.colNums()[d];

    if (lsa.master())
    {
        const label m = lsa.partSize(d);
        const label M = lsa.globalSize(d);
        const label nz = Foam::min(label(FullSType::nComponents),M);

        if (matrices_.set(d))
            MatDestroy(&matrices_[d]);

        matrices_.set(d, new Mat);

        Mat& mat = matrices_[d];

        MatCreate(comm, &mat);
        MatSetSizes(mat, m, m, M, M);

        if (Pstream::parRun() && lsa.nParts() > 1)
        {
            MatSetType(mat, MATMPIAIJ);
            MatMPIAIJSetPreallocation(mat, nz, NULL, nz, NULL);
        }
        else
        {
            MatSetType(mat, MATSEQAIJ);
            MatSeqAIJSetPreallocation(mat, nz, NULL);
        }

        MatSetUp(mat);

        int row = lsa.partStart(d);
        forAll(coeffs, proc)
        {
            forAll(coeffs[proc], i)
            {
                forAll(coeffs[proc][i], j)
                    if (colNums[proc][i][j] >= 0 && coeffs[proc][i][j] != 0)
                        MatSetValue
                        (
                            mat,
                            row,
                            colNums[proc][i][j],
                            coeffs[proc][i][j],
                            INSERT_VALUES
                        );

                row++;
            }
        }

        MatAssemblyBegin(mat, MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(mat, MAT_FINAL_ASSEMBLY);

        // View matrix for debugging:

        // MatView
        // (
        //     mat,
        //     PETSC_VIEWER_STDOUT_
        //     (
        //         Pstream::parRun()
        //       ? PstreamGlobals::lsaGetComm(lsa.masterCommNum())
        //       : PETSC_COMM_SELF
        //     )
        // );
    }
}

}

}

}
