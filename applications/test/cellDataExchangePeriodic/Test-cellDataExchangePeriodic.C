#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

#include "cellDataExchange.H"
#include "parallelBoundary.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class Type, class MeshType>
void testDataExchange(const fvMesh& fvMsh)
{
    meshField<Type,MeshType> f("f", fvMsh);

    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    const labelVector N(12,12,12);
    const labelVector start(fvMsh.msh().decomp().map().myBrickPartStart());

    for (int dir = 0; dir < 2; dir++)
    {
        forAllCells(f, d, i, j, k)
        {
            f(d,i,j,k) = cc(d,i,j,k)[dir]*pTraits<Type>::one;
        }

        f.correctBoundaryConditions();

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            // Add sample cells that are outside this processor domain

            DynamicList<labelVector> cells;

            forAll(fvMsh.msh().boundaries(), i)
            {
                const boundary& b = fvMsh.msh().boundaries()[i];

                if (b.castable<periodicBoundary>())
                {
                    const labelVector bo = b.offset();

                    const labelVector S(fvMsh.S<MeshType>(d,bo));
                    const labelVector E(fvMsh.E<MeshType>(d,bo));

                    labelVector ijk;

                    for (int j = 1; j <= 2; j++)
                    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                    {
                        cells.append(ijk + bo*j);
                    }
                }
            }

            // Also add some local cells

            cells.append(labelVector(1,1,0));
            cells.append(labelVector(2,0,0));
            cells.append(labelVector(0,2,1));
            cells.append(labelVector(1,2,2));

            cellDataExchange<MeshType> exchange(cells, fvMsh, 0, d);

            List<Type> data(move(exchange(f)));

            // Check the coordinate

            forAll(cells, i)
            {
                const labelVector ijk
                (
                    (start.x() + cells[i].x() + N.x()) % N.x(),
                    (start.y() + cells[i].y() + N.y()) % N.y(),
                    (start.z() + cells[i].z() + N.z()) % N.z()
                );

                const vector C
                (
                    cmptDivide
                    (
                        vector(ijk) + vector(unitXYZ)/2 + MeshType::shift[d],
                        vector(N)
                    )
                );


                if (Foam::mag(C[dir]*pTraits<Type>::one - data[i]) > 1e-8)
                    FatalErrorInFunction
                        << "Test 1 failed" << endl << abort(FatalError);
            }
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
