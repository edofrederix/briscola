#include "arguments.H"
#include "fv.H"
#include "immersedBoundary.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"

    // Test point 1 is outside of the IB
    vector testPoint1(0,0,0);
    // Test point 2 is inside of the IB
    vector testPoint2(0.5,0.5,0.5);
    // Test point 3 is inside cylinder2 but not cylinder1
    vector testPoint3(0.99,0.99,0.5);
    // Test point 3 is inside sphere1
    vector testPoint4(1,1,1);

    forAll(fvMsh.ibs<colocated>(), ib)
    {
        if (fvMsh.ibs<colocated>()[ib].isInside(testPoint1))
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

        if (!fvMsh.ibs<colocated>()[ib].isInside(testPoint2))
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

        if (!fvMsh.ibs<colocated>()[ib].isInside(testPoint3))
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

        if (!fvMsh.ibs<colocated>()[ib].isInside(testPoint4))
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


    forAll(fvMsh.ibs<colocated>(), ib)
    {
        if (fvMsh.ibs<colocated>()[ib].isInside(cc[0][0](testPoint5)))
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

        if (!fvMsh.ibs<colocated>()[ib].isInside(cc[0][0](testPoint6)))
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


    // Test wall distance functions

    // Line passes through sphere2 center
    vector c(10,10,11.12345);
    vector nb(10,10,10);

    forAll(fvMsh.ibs<colocated>(), ib)
    {
        if (mag(fvMsh.ibs<colocated>()[ib].wallDistance(c,nb) - 0.12345) > 1e-5)
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
    }

    // Line is tangent to the sphere2
    c = vector(10,11,11);
    nb = vector(10,10,11);

    forAll(fvMsh.ibs<colocated>(), ib)
    {
        if (mag(fvMsh.ibs<colocated>()[ib].wallDistance(c,nb) - 1.0) > 1e-5)
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
    }

    // Line arbitrarily crosses sphere2
    c = vector(10,12,10.5);
    nb = vector(10,10,10.5);

    forAll(fvMsh.ibs<colocated>(), ib)
    {
        if
        (
            mag
            (
                fvMsh.ibs<colocated>()[ib].wallDistance(c,nb)
              - (2.0 - Foam::sqrt(0.75))
            ) > 1e-5
        )
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
    }

    // Line passes through cylinder3 axis
    c = vector(11.12345,-1,0.5);
    nb = vector(10,-1,0.5);

    forAll(fvMsh.ibs<colocated>(), ib)
    {
        if (mag(fvMsh.ibs<colocated>()[ib].wallDistance(c,nb) - 0.12345) > 1e-5)
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
    }

    // Line is tangent to cylinder3
    c = vector(11,-2,0.5);
    nb = vector(10,-2,0.5);

    forAll(fvMsh.ibs<colocated>(), ib)
    {
        if (mag(fvMsh.ibs<colocated>()[ib].wallDistance(c,nb) - 1) > 1e-5)
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
    }

    // Line arbitrarily crosses cylinder3
    c = vector(12,-0.5,0.5);
    nb = vector(10,-0.5,0.5);

    forAll(fvMsh.ibs<colocated>(), ib)
    {
        if
        (
            mag
            (
                fvMsh.ibs<colocated>()[ib].wallDistance(c,nb)
              - (2.0 - Foam::sqrt(0.75))
            ) > 1e-5
        )
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
    }
}
