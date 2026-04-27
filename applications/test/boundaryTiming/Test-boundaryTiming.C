#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

#include "faceFields.H"
#include "stencilFields.H"
#include "diagStencilFields.H"

#include "ticToc.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

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

    initTicToc(32)

    colocatedScalarField coloScalar("coloScalar", fvMsh);
    colocatedVectorField coloVector("coloVector", fvMsh);

    staggeredScalarField stagScalar("stagScalar", fvMsh);
    staggeredVectorField stagVector("stagVector", fvMsh);

    coloScalar.makeDeep();
    coloVector.makeDeep();

    stagScalar.makeDeep();
    stagVector.makeDeep();

    forAllCells(coloScalar, l, d, i, j, k)
    {
        coloScalar(l,d,i,j,k) = scalar(l+d+i+j+k);
        coloVector(l,d,i,j,k) = scalar(l+d+i+j+k)*vector::one;
    }

    forAllCells(stagScalar, l, d, i, j, k)
    {
        stagScalar(l,d,i,j,k) = scalar(l+d+i+j+k);
        stagVector(l,d,i,j,k) = scalar(l+d+i+j+k)*vector::one;
    }

    for (label outer = 0; outer < 5; outer++)
    {
        // First level

        for (label iter = 0; iter < 5; iter++)
        {
            Tic(0)
            coloScalar[0].correctBoundaryConditions();
            Toc(0)

            Tic(1)
            coloVector[0].correctBoundaryConditions();
            Toc(1)

            Tic(2)
            stagScalar[0].correctBoundaryConditions();
            Toc(2)

            Tic(3)
            stagVector[0].correctBoundaryConditions();
            Toc(3)
        }

        // Second level

        for (label iter = 0; iter < 5; iter++)
        {
            Tic(4)
            coloScalar[1].correctBoundaryConditions();
            Toc(4)

            Tic(5)
            coloVector[1].correctBoundaryConditions();
            Toc(5)

            Tic(6)
            stagScalar[1].correctBoundaryConditions();
            Toc(6)

            Tic(7)
            stagVector[1].correctBoundaryConditions();
            Toc(7)
        }

        // Third level

        for (label iter = 0; iter < 5; iter++)
        {
            Tic(8)
            coloScalar[2].correctBoundaryConditions();
            Toc(8)

            Tic(9)
            coloVector[2].correctBoundaryConditions();
            Toc(9)

            Tic(10)
            stagScalar[2].correctBoundaryConditions();
            Toc(10)

            Tic(11)
            stagVector[2].correctBoundaryConditions();
            Toc(11)
        }

        // Whole field

        for (label iter = 0; iter < 5; iter++)
        {
            Tic(12)
            coloScalar.correctBoundaryConditions();
            Toc(12)

            Tic(13)
            coloVector.correctBoundaryConditions();
            Toc(13)

            Tic(14)
            stagScalar.correctBoundaryConditions();
            Toc(14)

            Tic(15)
            stagVector.correctBoundaryConditions();
            Toc(15)
        }
    }

    Info<< "First level" << nl
        << "    Colocated scalar = " << ticTocs[0]/1e3 << " ms" << nl
        << "    Colocated vector = " << ticTocs[1]/1e3 << " ms" << nl
        << "    Staggered scalar = " << ticTocs[2]/1e3 << " ms" << nl
        << "    Staggered vector = " << ticTocs[3]/1e3 << " ms" << endl;

    Info<< "Second level" << nl
        << "    Colocated scalar = " << ticTocs[4]/1e3 << " ms" << nl
        << "    Colocated vector = " << ticTocs[5]/1e3 << " ms" << nl
        << "    Staggered scalar = " << ticTocs[6]/1e3 << " ms" << nl
        << "    Staggered vector = " << ticTocs[7]/1e3 << " ms" << endl;

    Info<< "Third level" << nl
        << "    Colocated scalar = " << ticTocs[8]/1e3 << " ms" << nl
        << "    Colocated vector = " << ticTocs[9]/1e3 << " ms" << nl
        << "    Staggered scalar = " << ticTocs[10]/1e3 << " ms" << nl
        << "    Staggered vector = " << ticTocs[11]/1e3 << " ms" << endl;

    Info<< "Whole field" << nl
        << "    Colocated scalar = " << ticTocs[12]/1e3 << " ms" << nl
        << "    Colocated vector = " << ticTocs[13]/1e3 << " ms" << nl
        << "    Staggered scalar = " << ticTocs[14]/1e3 << " ms" << nl
        << "    Staggered vector = " << ticTocs[15]/1e3 << " ms" << endl;

    Info<< endl;
}
