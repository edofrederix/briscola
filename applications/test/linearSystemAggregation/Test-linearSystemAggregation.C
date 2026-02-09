#include "arguments.H"
#include "Time.H"

#include "linearSystemAggregation.H"
#include "linearSystem.H"
#include "imSchemes.H"
#include "solver.H"
#include "EigenLinearSystem.H"

#include "OFstream.H"
#include "IFstream.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

template<class SType, class Type, class MeshType>
void test(const fvMesh& fvMsh)
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

    restrict(sys.x());
    restrict(sys.b());

    for (int nParts = 1; nParts <= Pstream::nProcs(); nParts++)
    forAll(fvMsh.msh(), l)
    if (nParts <= fvMsh[l].decomp().members().size())
    {
        linearSystemAggregation<SType,Type,MeshType> lsa
        (
            sys,
            l,
            nParts
        );

        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            const List<List<FixedList<label,SType::nCsComponents>>>&
                colNums = lsa.colNums()[d];

            List<List<SType>> rows;
            lsa.rowCoeffs(rows, sys, d);

            List<Type> rhs;
            lsa.rhsSource(rhs, sys, d);

            scalarList values;
            labelList inners;
            labelList outers;
            lsa.compressedRowFormat(values, inners, outers, sys, d);
            lsa.compressedRowFormat(values, inners, outers, sys, d, true);

            // Collect global indices and coefficients

            List<DynamicList<label>> gRows(Pstream::nProcs());
            List<DynamicList<label>> gCols(Pstream::nProcs());
            List<DynamicList<scalar>> coeffs(Pstream::nProcs());

            const label start = lsa.starts()[d];

            label cur = 0;
            forAll(rows, proc)
            {
                forAll(rows[proc], row)
                {
                    const label gRow = start + cur;

                    for (label col = 0; col < SType::nCsComponents; col++)
                    {
                        const label gCol = colNums[proc][row][col];

                        if (gCol > -1)
                        {
                            const scalar coeff = rows[proc][row][col];

                            gRows[Pstream::myProcNo()].append(gRow);
                            gCols[Pstream::myProcNo()].append(gCol);
                            coeffs[Pstream::myProcNo()].append(coeff);
                        }
                    }

                    cur++;
                }
            }

            Pstream::gatherList(gRows);
            Pstream::gatherList(gCols);
            Pstream::gatherList(coeffs);

            // Write to a file, so that we can check if linear systems
            // remain independent of the number of parts

            if (Pstream::master())
            {
                const word fileName
                (
                    f.name() + "_" +
                    SType::typeName + "_" +
                    Foam::name(nParts) + "_" +
                    Foam::name(l) + "_" +
                    Foam::name(d)
                );

                OFstream os(fileName);

                DynamicList<label> gRowsLin;
                DynamicList<label> gColsLin;
                DynamicList<scalar> coeffsLin;

                forAll(gRows, proc)
                {
                    forAll(gRows[proc], i)
                    {
                        gRowsLin.append(gRows[proc][i]);
                        gColsLin.append(gCols[proc][i]);
                        coeffsLin.append(coeffs[proc][i]);
                    }
                }

                os << gRowsLin << nl;
                os << gColsLin << nl;
                os << coeffsLin << nl;
            }
        }

        // Test Eigen linear system which uses LSA

        EigenLinearSystem<SType,Type,MeshType> E
        (
            sys,
            l,
            nParts
        );
    }

    // Check if systems are independent of the number of parts

    if (Pstream::master())
    forAll(fvMsh.msh(), l)
    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        labelList gRows(Pstream::nProcs());
        labelList gCols(Pstream::nProcs());
        scalarList coeffs(Pstream::nProcs());

        for (int nParts = 1; nParts <= Pstream::nProcs(); nParts++)
        {
            if (nParts <= fvMsh[l].decomp().members().size())
            {
                const word fileName
                (
                    f.name() + "_" +
                    SType::typeName + "_" +
                    Foam::name(nParts) + "_" +
                    Foam::name(l) + "_" +
                    Foam::name(d)
                );

                IFstream is(fileName);

                if (nParts == 1)
                {
                    is >> gRows;
                    is >> gCols;
                    is >> coeffs;
                }
                else
                {
                    labelList gRowsNext(Pstream::nProcs());
                    labelList gColsNext(Pstream::nProcs());
                    scalarList coeffsNext(Pstream::nProcs());

                    is >> gRowsNext;
                    is >> gColsNext;
                    is >> coeffsNext;

                    forAll(gRowsNext, i)
                        if (gRowsNext[i] != gRows[i])
                            FatalErrorInFunction
                                << "Error 1" << abort(FatalError);

                    forAll(gColsNext, i)
                        if (gColsNext[i] != gCols[i])
                            FatalErrorInFunction
                                << "Error 2" << abort(FatalError);

                    forAll(coeffsNext, i)
                        if (coeffsNext[i] != coeffsNext[i])
                            FatalErrorInFunction
                                << "Error 3" << abort(FatalError);
                }
            }
        }
    }
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    test<stencil,scalar,colocated>(fvMsh);
    test<stencil,vector,colocated>(fvMsh);

    test<stencil,scalar,staggered>(fvMsh);
    test<stencil,vector,staggered>(fvMsh);
}
