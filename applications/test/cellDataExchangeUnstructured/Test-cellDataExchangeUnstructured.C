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

inline scalar coord(const labelVector& ijk)
{
    return (ijk.z()+0.5)/16.0;
}

template<class Type>
inline bool test(const Type& value, const labelVector& ijk)
{
    return
        Foam::mag
        (
            value
          - coord(ijk)*pTraits<Type>::one
        )
      < 1e-8;
}

template<class Type>
void testDataExchange(const fvMesh& fvMsh)
{
    meshField<Type,colocated> f("f", fvMsh);

    const meshField<vector,colocated>& cc =
        fvMsh.metrics<colocated>().cellCenters();

    forAllDirections(f, d, i, j, k)
    {
        f(d,i,j,k) = cc(d,i,j,k).z()*pTraits<Type>::one;
    }

    f.correctBoundaryConditions();

    List<labelVector> cells;

    forAll(f.boundaryConditions(), i)
    {
        const boundaryCondition<Type,colocated>& bc =
            f.boundaryConditions()[i];

        if (bc.baseType() == PARALLELBC)
        {
            const labelVector S(bc.S(0,0));
            const labelVector E(bc.E(0,0));

            const labelVector bo = bc.boundaryOffset();

            labelVector ijk;

            // Add some cells across this parallel boundary

            int count = 0;
            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()+=3)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()+=4)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()+=5)
            {
                cells.append(ijk + bo*(1+1*((count++)%3)));
            }
        }
    }

    cellDataExchange<colocated> exchange(cells, fvMsh, 0);

    List<Type> data(exchange(f));

    // Check the coordinate

    forAll(data, i)
        if (!test<Type>(data[i],cells[i]))
            FatalErrorInFunction
                << "Test 1a failed" << endl << abort(FatalError);
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

    testDataExchange<scalar>(fvMsh);
    testDataExchange<vector>(fvMsh);
    testDataExchange<tensor>(fvMsh);
    testDataExchange<sphericalTensor>(fvMsh);
    testDataExchange<symmTensor>(fvMsh);
    testDataExchange<diagTensor>(fvMsh);
}
