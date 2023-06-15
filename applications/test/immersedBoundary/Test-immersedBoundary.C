#include "arguments.H"
#include "fv.H"
#include "immersedBoundary.H"
#include "cylinder.H"

using namespace Foam;
using namespace briscola;
using namespace fv;
using namespace ibm;


int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    IOdictionary solverDict
    (
        IOobject
        (
            "briscolaStaggeredDict",
            fvMsh.time().system(),
            fvMsh.time(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    immersedBoundary IB(solverDict, fvMsh);

    // Test point 1 is outside of the IB
    vector testPoint1(0,0,0);
    // Test point 2 is inside of the IB
    vector testPoint2(0.5,0.5,0.5);
    // Test point 3 is inside cylinder2 but not cylinder1
    vector testPoint3(0.99,0.99,0.5);
    // Test point 3 is inside sphere1
    vector testPoint4(1,1,1);

    if (IB.isInside(testPoint1))
    {
        FatalError
            << "IB test 1 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB test 1 successful!" << endl;
    }

    if (!IB.isInside(testPoint2))
    {
        FatalError
            << "IB test 2 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB test 2 successful!" << endl;
    }

    if (!IB.isInside(testPoint3))
    {
        FatalError
            << "IB test 3 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB test 3 successful!" << endl;
    }

    if (!IB.isInside(testPoint4))
    {
        FatalError
            << "IB test 4 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB test 4 successful!" << endl;
    }



    // Test point 5 is defined by its cell indices
    // and is outside the IB
    labelVector testPoint5(1,1,1);

    // Test point 5 is defined by its cell indices
    // and is inside the IB
    labelVector testPoint6(64,64,64);

    // Cell centers
    // cc[0][0] returns the collocated cell center
    // locations on the highest mesh level
    const staggeredVectorField& cc =
        fvMsh.metrics<staggered>().cellCenters();


    if (IB.isInside(cc[0][0](testPoint5)))
    {
        FatalError
            << "IB test 5 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB test 5 successful!" << endl;
    }

    if (!IB.isInside(cc[0][0](testPoint6)))
    {
        FatalError
            << "IB test 6 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB test 6 successful!" << endl;
    }

}
