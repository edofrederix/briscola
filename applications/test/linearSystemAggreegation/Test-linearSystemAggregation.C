#include "arguments.H"
#include "Time.H"

#include "linearSystemAggregation.H"
#include "linearSystem.H"
#include "imSchemes.H"
#include "solver.H"
#include "EigenLinearSystem.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class SType, class Type, class MeshType>
void test(const fvMesh& fvMsh, const label nParts)
{
    meshField<Type,MeshType> f
    (
        "f-"
      + word(MeshType::typeName)
      + "-"
      + word(pTraits<Type>::typeName),
        fvMsh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE,
        true
    );

    f = 0.1*pTraits<Type>::one;

    linearSystem<SType,Type,MeshType> sys(im::laplacian<SType>(f));
    sys -= im::ddt(f);

    sys.singular();
    sys.diagonal();

    sys.x().restrict();
    sys.b().restrict();

    forAll(fvMsh.msh(), l)
    {
        for (int i = 0; i < 2; i++)
        {
            const List<bool> singular(MeshType::numberOfDirections, Switch(i));

            linearSystemAggregation<SType,Type,MeshType> lsa
            (
                sys,
                l,
                nParts
            );

            for (int d = 0; d < MeshType::numberOfDirections; d++)
            {
                List<List<typename SType::fullStencilType>> rows;
                lsa.rowCoeffs(rows, sys, d);

                List<Type> rhs;
                lsa.rhsSource(rhs, sys, d);

                scalarList values;
                labelList inners;
                labelList outers;
                lsa.compressedRowFormat(values, inners, outers, sys, d);
                lsa.compressedRowFormat(values, inners, outers, sys, d, true);
            }

            // Test Eigen linear system uses LSA

            EigenLinearSystem<SType,Type,MeshType> E
            (
                sys,
                l,
                nParts
            );
        }
    }
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    for (int nParts = 1; nParts <= Pstream::nProcs(); nParts++)
    {
        test<stencil,scalar,colocated>(fvMsh, nParts);
        test<symmStencil,scalar,colocated>(fvMsh, nParts);

        test<stencil,vector,colocated>(fvMsh, nParts);
        test<symmStencil,vector,colocated>(fvMsh, nParts);

        test<stencil,scalar,staggered>(fvMsh, nParts);
        test<symmStencil,scalar,staggered>(fvMsh, nParts);

        test<stencil,vector,staggered>(fvMsh, nParts);
        test<symmStencil,vector,staggered>(fvMsh, nParts);
    }
}
