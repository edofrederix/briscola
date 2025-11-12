#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

#include "pointDataExchange.H"

#include "randomGenerator.H"
#include "constants.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int seed = 0;
const int N = 100;

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
    const meshField<vector,MeshType>& cc =
        fvMsh.metrics<MeshType>().cellCenters();

    meshField<Type,MeshType> f("f", fvMsh);
    f.makeDeep();
    f = Zero;

    const faceScalar bb(fvMsh.msh().boundingBox());

    for (int dir = 0; dir < 3; dir++)
    {
        // Copy all data including ghosts, not just cell data

        forAll(f, l)
            forAll(f[l], d)
                forAllBlock(f[l][d].B(), i, j, k)
                    f[l][d].B()(i,j,k) =
                        cc[l][d].B()(i,j,k)[dir]*pTraits<Type>::one;

        forAll(f, l)
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            // Generate random points, assuming a domain that coincides with its
            // bounding box

            randomGenerator R(Pstream::myProcNo() + seed++);

            vectorList points(N);
            forAll(points, c)
                points[c] = vector
                (
                    R.scalarAB(bb.left(),   bb.right()),
                    R.scalarAB(bb.bottom(), bb.top()  ),
                    R.scalarAB(bb.aft(),    bb.fore() )
                );

            pointDataExchange<MeshType> exchange(points, fvMsh, l, d);

            List<Type> data(move(exchange(f)));

            // Check the coordinate

            forAll(data, i)
                if
                (
                    Foam::mag
                    (
                        Foam::mag(data[i])/Foam::mag(pTraits<Type>::one)
                      - points[i][dir]
                    ) > 1e-3
                )
                    FatalErrorInFunction
                        << "Test 1 failed" << endl << abort(FatalError);
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
