#include "arguments.H"
#include "Time.H"
#include "fvMesh.H"
#include "pointInterpolator.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class MeshType>
void testCellCenters(const fvMesh& fvMsh)
{
    const meshField<vector,MeshType>& cc =
        fvMsh.template metrics<MeshType>().cellCenters();

    forAllCells(cc, l, d, i, j, k)
    {
        // Near edges, move the point a little bit inward to avoid issues

        const labelVector shift
        (
            i == 0 ? 1 : i == cc.I(l,d).right()-1 ? -1 : 0,
            j == 0 ? 1 : j == cc.I(l,d).top()  -1 ? -1 : 0,
            k == 0 ? 1 : k == cc.I(l,d).fore() -1 ? -1 : 0
        );

        const vector point =
            (1.0-1e-12)*cc(l,d,i,j,k)
          + 1e-12*cc(l,d,i+shift.x(),j+shift.y(),k+shift.z());

        labelVector ijk = fvMsh.findCell<MeshType>(point, l, d);

        if (ijk != labelVector(i,j,k))
            FatalErrorInFunction
                << "Test 1 failed" << endl << abort(FatalError);
    }
}

template<class Type, class MeshType>
void testInterpolations(const fvMesh& fvMsh)
{
    const meshField<vector,MeshType>& cc =
        fvMsh.template metrics<MeshType>().cellCenters();

    meshField<Type,MeshType> field("field", fvMsh);
    field.makeDeep();

    field = Zero;

    forAllCells(field, l, d, i, j, k)
    {
        field(l,d,i,j,k) = (i+j+k+l+d)*pTraits<Type>::one;
    }

    forAll(cc, l)
    {
        forAll(cc[l], d)
        {
            const meshDirection<vector,MeshType>& ccld = cc[l][d];

            vectorList points(ccld.size());

            label c = 0;
            forAllCells(ccld, i, j, k)
            {
                // Near edges, move the point a little bit inward to avoid
                // issues

                const labelVector shift
                (
                    i == 0 ? 1 : i == cc.I(l,d).right()-1 ? -1 : 0,
                    j == 0 ? 1 : j == cc.I(l,d).top()  -1 ? -1 : 0,
                    k == 0 ? 1 : k == cc.I(l,d).fore() -1 ? -1 : 0
                );

                points[c++] =
                    (1.0-1e-12)*ccld(i,j,k)
                  + 1e-12*ccld(i+shift.x(),j+shift.y(),k+shift.z());
            }

            autoPtr<pointInterpolator<MeshType>> nInterpPtr
            (
                pointInterpolator<MeshType>::New
                (
                    fvMsh,
                    points,
                    "nearest",
                    false,
                    l,
                    d
                )
            );

            autoPtr<pointInterpolator<MeshType>> lInterpPtr
            (
                pointInterpolator<MeshType>::New
                (
                    fvMsh,
                    points,
                    "linear",
                    false,
                    l,
                    d
                )
            );

            pointInterpolator<MeshType>& nInterp = nInterpPtr();
            pointInterpolator<MeshType>& lInterp = lInterpPtr();

            List<Type> nValues(nInterp(field));
            List<Type> lValues(lInterp(field));

            c = 0;
            forAllCells(ccld, i, j, k)
                if (Foam::mag(nValues[c++] - field(l,d,i,j,k)) > 1e-4)
                    FatalErrorInFunction
                        << "test 2 failed" << endl << abort(FatalError);

            c = 0;
            forAllCells(ccld, i, j, k)
                if (Foam::mag(lValues[c++] - field(l,d,i,j,k)) > 1e-4)
                    FatalErrorInFunction
                        << "test 3 failed" << endl << abort(FatalError);
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

    // Find cell centers and see if they match

    testCellCenters<colocated>(fvMsh);

    if (fvMsh.structured())
        testCellCenters<staggered>(fvMsh);

    // Interpolate points

    testInterpolations<scalar,colocated>(fvMsh);
    testInterpolations<vector,colocated>(fvMsh);
    testInterpolations<tensor,colocated>(fvMsh);
    testInterpolations<sphericalTensor,colocated>(fvMsh);
    testInterpolations<symmTensor,colocated>(fvMsh);
    testInterpolations<diagTensor,colocated>(fvMsh);

    if (fvMsh.structured())
    {
        testInterpolations<scalar,staggered>(fvMsh);
        testInterpolations<vector,staggered>(fvMsh);
        testInterpolations<tensor,staggered>(fvMsh);
        testInterpolations<sphericalTensor,staggered>(fvMsh);
        testInterpolations<symmTensor,staggered>(fvMsh);
        testInterpolations<diagTensor,staggered>(fvMsh);
    }
}
