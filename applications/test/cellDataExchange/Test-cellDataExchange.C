#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

#include "cellDataExchange.H"
#include "parallelPartPatch.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class MeshType>
inline scalar coord
(
    const label d,
    const labelVector& ijk,
    const label coordDir
)
{
    const vector c =
        vector(MeshType::shift[d])/64.0
      + 0.5*vector(unitX)*(Pstream::myProcNo() > 1)
      + 0.5*vector(unitY)*(Pstream::myProcNo() % 2)
      + (vector(ijk)+0.5*vector(unitXYZ))/64.0;

    return c[coordDir];
}

template<class Type, class MeshType>
inline bool test
(
    const Type& value,
    const labelVector& ijk,
    const label d,
    const label coordDir
)
{
    return
        Foam::mag
        (
            value
          - coord<MeshType>(d,ijk,coordDir)*pTraits<Type>::one
        )
      < 1e-8;
}

template<class Type, class MeshType>
void testDataExchange(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> f("f", fvMsh);

    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    for (int dir = 0; dir < 2; dir++)
    {
        forAllDirections(f, d, i, j, k)
        {
            f(d,i,j,k) = cc(d,i,j,k)[dir]*pTraits<Type>::one;
        }

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const faceLabel I = fvMsh.I<MeshType>(d);

            // Add some sample cells that are outside this processor domain

            List<labelVector> cells;

            forAll(fvMsh.msh().partPatches(), i)
            {
                const partPatch& patch = fvMsh.msh().partPatches()[i];

                if (patch.typeNum() == parallelPartPatch::typeNumber)
                {
                    const labelVector bo = patch.boundaryOffset();

                    for (int j = 0; j < 4+Pstream::myProcNo(); j++)
                        cells.append
                        (
                            (
                                bo.x() == 0 ? unitX*(j*3+2)
                              : bo.x() == 1 ? (I.right()-1)*unitX
                              :               zeroXYZ
                            )
                          + (
                                bo.y() == 0 ? unitY*(j*3+2)
                              : bo.y() == 1 ? (I.top()-1)*unitY
                              :               zeroXYZ
                            )
                          + labelVector(bo*(j+1))
                        );
                }
            }

            cellDataExchange<MeshType> exchange(cells, fvMsh, d);

            List<Type> data(move(exchange(f)));

            // Check the coordinate

            forAll(data, i)
                if (!test<Type,MeshType>(data[i],cells[i],d,dir))
                    FatalErrorInFunction
                        << "Test 1a failed" << endl << abort(FatalError);
        }
    }
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"

    IOdictionary meshDict
    (
        IOobject
        (
            runTime.system()/"briscolaMeshDict",
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    fvMesh fvMsh(meshDict, runTime);

    testDataExchange<scalar,colocated>(fvMsh);
    testDataExchange<vector,colocated>(fvMsh);
    testDataExchange<tensor,colocated>(fvMsh);
    testDataExchange<sphericalTensor,colocated>(fvMsh);
    testDataExchange<symmTensor,colocated>(fvMsh);
    testDataExchange<diagTensor,colocated>(fvMsh);

    testDataExchange<scalar,staggered>(fvMsh);
    testDataExchange<vector,staggered>(fvMsh);
    testDataExchange<tensor,staggered>(fvMsh);
    testDataExchange<sphericalTensor,staggered>(fvMsh);
    testDataExchange<symmTensor,staggered>(fvMsh);
    testDataExchange<diagTensor,staggered>(fvMsh);
}
