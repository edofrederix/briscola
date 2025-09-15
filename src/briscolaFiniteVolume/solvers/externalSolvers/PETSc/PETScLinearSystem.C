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
    matrices_(MeshType::numberOfDirections, nullptr),
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
    const linearSystemAggregation<SType,Type,MeshType>& lsa = *this;

    const MPI_Comm& comm =
        Pstream::parRun()
      ? PstreamGlobals::lsaGetComm(lsa.masterCommNum())
      : PETSC_COMM_SELF;

    List<List<SType>> coeffs;
    lsa.rowCoeffs(coeffs, sys_, d);

    List<List<FixedList<label,SType::nComponents>>> colNums =
        lsa.colNums()[d];

    if (lsa.master())
    {
        const label m = lsa.partSize(d);
        const label M = lsa.globalSize(d);
        const label nz = Foam::min(label(SType::nComponents),M);

        if (matrices_[d])
        {
            MatDestroy(&matrices_[d]);
            matrices_[d] = nullptr;
        }

        MatCreate(comm, &matrices_[d]);
        MatSetSizes(matrices_[d], m, m, M, M);

        if (Pstream::parRun() && lsa.nParts() > 1)
        {
            MatSetType(matrices_[d], MATMPIAIJ);
            MatMPIAIJSetPreallocation(matrices_[d], nz, NULL, nz, NULL);
        }
        else
        {
            MatSetType(matrices_[d], MATSEQAIJ);
            MatSeqAIJSetPreallocation(matrices_[d], nz, NULL);
        }

        MatSetUp(matrices_[d]);

        int row = lsa.partStart(d);
        forAll(coeffs, proc)
        {
            forAll(coeffs[proc], i)
            {
                forAll(coeffs[proc][i], j)
                    if (colNums[proc][i][j] >= 0 && coeffs[proc][i][j] != 0)
                        MatSetValue
                        (
                            matrices_[d],
                            row,
                            colNums[proc][i][j],
                            coeffs[proc][i][j],
                            INSERT_VALUES
                        );

                row++;
            }
        }

        MatAssemblyBegin(matrices_[d], MAT_FINAL_ASSEMBLY);
        MatAssemblyEnd(matrices_[d], MAT_FINAL_ASSEMBLY);

        // View matrix for debugging:

        // MatView
        // (
        //     matrices_[d],
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
