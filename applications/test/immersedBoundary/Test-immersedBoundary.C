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

    immersedBoundary<colocated> IB(solverDict, fvMsh);

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
    // cc[0][0] returns the colocated cell center
    // locations on the highest mesh level
    const colocatedVectorField& cc =
        fvMsh.metrics<colocated>().cellCenters();


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


    // Test wall distance functions

    sphere sphere2(vector(10,10,10), 1, false);

    // Line passes through sphere center
    vector c(10,10,11.12345);
    vector nb(10,10,10);

    if (mag(sphere2.wallDistance(c,nb) - 0.12345) > 1e-5)
    {
        FatalError
            << "IB wall distance test 1 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB wall distance test 1 successful!" << endl;
    }

    // Line is tangent to the sphere
    c = vector(10,11,11);
    nb = vector(10,10,11);

    if (mag(sphere2.wallDistance(c,nb) - 1.0) > 1e-5)
    {
        FatalError
            << "IB wall distance test 2 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB wall distance test 2 successful!" << endl;
    }

    // Line arbitrarily crosses sphere
    c = vector(10,12,10.5);
    nb = vector(10,10,10.5);

    if (mag(sphere2.wallDistance(c,nb) - (2.0 - Foam::sqrt(0.75))) > 1e-5)
    {
        FatalError
            << "IB wall distance test 3 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB wall distance test 3 successful!" << endl;
    }

    cylinder cylinder3(vector(0,0,0), vector(0,0,1), 1.0, false);

    // Line passes through cylinder axis
    c = vector(1.12345,0,0.5);
    nb = vector(0,0,0.5);

    if (mag(cylinder3.wallDistance(c,nb) - 0.12345) > 1e-5)
    {
        FatalError
            << "IB wall distance test 4 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB wall distance test 4 successful!" << endl;
    }

    // Line is tangent to the cylinder
    c = vector(1,1,0.5);
    nb = vector(0,1,0.5);

    if (mag(cylinder3.wallDistance(c,nb) - 1) > 1e-5)
    {
        FatalError
            << "IB wall distance test 5 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB wall distance test 5 successful!" << endl;
    }

    // Line arbitrarily crosses cylinder
    c = vector(2,0.5,0.5);
    nb = vector(0,0.5,0.5);

    if (mag(cylinder3.wallDistance(c,nb) - (2.0 - Foam::sqrt(0.75))) > 1e-5)
    {
        FatalError
            << "IB wall distance test 6 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB wall distance test 6 successful!" << endl;
    }

    // Line is parallel to cylinder axis and intersects end caps
    c = vector(0.5,0.5,2);
    nb = vector(0.5,0.5,0.5);

    if (mag(cylinder3.wallDistance(c,nb) - 1) > 1e-5)
    {
        FatalError
            << "IB wall distance test 7 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB wall distance test 7 successful!" << endl;
    }

    // Line intersects cylinder end cap
    c = vector(0,1,2);
    nb = vector(0,0,1);

    if (mag(cylinder3.wallDistance(c,nb) - Foam::sqrt(2.0)) > 1e-5)
    {
        FatalError
            << "IB wall distance test 8 failed."
            << endl;
        FatalError.exit();
    }
    else
    {
        Info << "IB wall distance test 8 successful!" << endl;
    }
}