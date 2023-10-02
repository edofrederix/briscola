#include "arguments.H"
#include "Time.H"
#include "fvMesh.H"
#include "pointInterpolator.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class Type>
void testColocatedInterpolations(const fvMesh& fvMsh)
{
    const partLevel& lvl = fvMsh.msh()[0];

    const meshDirection<vector,colocated>& cc =
        fvMsh.template metrics<colocated>().cellCenters().direction();

    vectorList points(lvl.size());

    label c = 0;
    forAllBlock(lvl, i, j, k)
    {
        points[c++] = cc(i,j,k);
    }

    autoPtr<pointInterpolator> nInterpPtr
    (
        pointInterpolator::New(fvMsh,points,"nearest")
    );

    autoPtr<pointInterpolator> lInterpPtr
    (
        pointInterpolator::New(fvMsh,points,"nearest")
    );

    pointInterpolator& nInterp = nInterpPtr();
    pointInterpolator& lInterp = lInterpPtr();

    meshField<Type,colocated> field("field", fvMsh);
    field = Zero;

    forAllLevels(field, l, d, i, j, k)
    {
        field(l,d,i,j,k) = (i+j+k+l+d)*pTraits<Type>::one;
    }

    List<Type> nValues(nInterp(field));
    List<Type> lValues(lInterp(field));

    c = 0;
    forAllBlock(lvl, i, j, k)
    {
        if (nValues[c] != field(i,j,k))
            FatalErrorInFunction
                << "test 3 failed" << endl << abort(FatalError);

        c++;
    }

    c = 0;
    forAllBlock(lvl, i, j, k)
    {
        if (lValues[c] != field(i,j,k))
            FatalErrorInFunction
                << "test 4 failed" << endl << abort(FatalError);

        c++;
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

    const mesh& msh = fvMsh.msh();

    // Find cell centers and see if they match

    forAll(msh, l)
    {
        const partLevel& lvl = msh[l];

        forAllBlock(lvl, i, j, k)
        {
            const vertexVector vertices(lvl.points().cellPoints(i,j,k));

            vector point = Zero;
            for(int v = 0; v < 8; v++)
                point += vertices[v]/8.0;

            labelVector ijk = msh.findCell(point, l);

            if (ijk != labelVector(i,j,k))
            {
                FatalErrorInFunction
                    << "Test 1 failed" << endl << abort(FatalError);
            }

            if (!lvl.points().pointInCell(point, ijk))
                FatalErrorInFunction
                    << "Test 2 failed" << endl << abort(FatalError);
        }
    }

    // Interpolate points

    testColocatedInterpolations<scalar>(fvMsh);
    testColocatedInterpolations<vector>(fvMsh);
    testColocatedInterpolations<tensor>(fvMsh);
    testColocatedInterpolations<sphericalTensor>(fvMsh);
    testColocatedInterpolations<symmTensor>(fvMsh);
    testColocatedInterpolations<diagTensor>(fvMsh);
    testColocatedInterpolations<faceScalar>(fvMsh);
    testColocatedInterpolations<faceVector>(fvMsh);
}
