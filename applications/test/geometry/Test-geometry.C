#include "arguments.H"

#include "geometry.H"
#include "constants.H"

using namespace Foam;
using namespace briscola;

const label Ntheta = 9;

template<class Type>
void test(const vertexVector v)
{
    // Test constructors

    Type h1(v);
    Type h2
    (
        v.lba(),
        v.rba(),
        v.lta(),
        v.rta(),
        v.lbf(),
        v.rbf(),
        v.ltf(),
        v.rtf()
    );

    // Calculate volumes

    if (h1.volume() < 1e-12)
        FatalErrorInFunction
            << "Invalid volume" << nl << abort(FatalError);

    if (tessellation(h1).volume() < 1e-12)
        FatalErrorInFunction
            << "Invalid volume" << nl << abort(FatalError);

    // Overlap of decomposition tets should be zero

    const tessellation T(h1);

    forAll(T, i)
    forAll(T, j)
    if (i != j)
    {
        const tetra& t0 = T[i];
        const tetra& t1 = T[j];

        if (t0.intersect(t1).size() > 0)
            FatalErrorInFunction
                << "Invalid tessellation" << nl << abort(FatalError);

        if (t0.intersect(t1).volume() > 1e-12)
            FatalErrorInFunction
                << "Invalid tessellation" << nl << abort(FatalError);
    }

    // Overlap of decomposition tets with self should be the same

    forAll(T, i)
    {
        const tetra& t = T[i];

        if ((t.intersect(t).volume() -t.volume()) > 1e-12)
            FatalErrorInFunction
                << "Invalid tessellation" << nl << abort(FatalError);
    }

    // Compare center of type and center of tessellation of type

    piped p1
    (
        h1.lba(),
        h1.rba() - h1.lba(),
        h1.lta() - h1.lba(),
        h1.lbf() - h1.lba()
    );

    const bool isPiped =
        Foam::mag(hexa(p1).v() - hexa(h1).v()) < 1e-12*vertexScalar::one;

    if (isPiped)
    {
        // Only of the object is a piped we can expect exact agreement

        if (Foam::mag(h1.center() - tessellation(h1).center()) > 1e-12)
            FatalErrorInFunction
                << "Invalid centers" << nl << abort(FatalError);
    }
    else
    {
        // Otherwise, let's check that the centers are at least within 10% of
        // each other

        if
        (
            Foam::mag(h1.center() - tessellation(h1).center())
          > Foam::cbrt(h1.volume())*0.1
        )
            FatalErrorInFunction
                << "Invalid centers" << nl << abort(FatalError);
    }

    // Normals

    vectorList ns(8);

    ns[0] = vector(1, 0, 0);
    ns[1] = vector(0, 1, 0);
    ns[2] = vector(0, 0, 1);
    ns[3] = -ns[0];
    ns[4] = -ns[1];
    ns[5] = -ns[2];
    ns[6] = vector(4.0/5.0, 1.0/5.0, 2.0*Foam::sqrt(2.0)/5.0);
    ns[7] = -ns[6];

    // Truncate

    forAll(ns, i)
    {
        const vector n1 = ns[i];
        const vector n2 = -ns[i];

        const vector c = h1.center();

        for (label j = -5; j <= 5; j++)
        {
            // Shifted center
            const vector d = c + n1*scalar(j)/10.0;

            const scalar C1 = - (d & n1);
            const scalar C2 = - (d & n2);

            const tessellation T1(h1.truncate(n1,C1));
            const tessellation T2(h1.truncate(n2,C2));

            tessellation T(T1);
            T.append(T2);

            // Two opposite truncations should sum into the same original volume

            if (Foam::mag(T1.volume() + T2.volume() - h1.volume()) > 1e-12)
                FatalErrorInFunction
                    << "Invalid truncations" << nl << abort(FatalError);

            if (Foam::mag(T.volume() - h1.volume()) > 1e-12)
                FatalErrorInFunction
                    << "Invalid truncations" << nl << abort(FatalError);
        }
    }
}

template<class Type>
void testRotated(const vertexVector v)
{
    const scalar theta =
        2.0*Foam::constant::mathematical::pi/scalar(Ntheta);

    // Rotate around x

    tensor T
    (
        vector(1,0,0),
        vector(0,Foam::cos(theta),-Foam::sin(theta)),
        vector(0,Foam::sin(theta),Foam::cos(theta))
    );

    vertexVector w = v;

    for (int i = 0; i < Ntheta; i++)
    {
        test<Type>(w);

        w = (T & w);
    }

    // Rotate around y

    T = tensor
    (
        vector(Foam::cos(theta),0,Foam::sin(theta)),
        vector(0,1,0),
        vector(-Foam::sin(theta),0,Foam::cos(theta))
    );

    w = v;

    for (int i = 0; i < Ntheta; i++)
    {
        test<Type>(w);

        w = (T & w);
    }

    // Rotate around z

    T = tensor
    (
        vector(Foam::cos(theta),-Foam::sin(theta),0),
        vector(Foam::sin(theta),Foam::cos(theta),0),
        vector(0,0,1)
    );

    w = v;

    for (int i = 0; i < Ntheta; i++)
    {
        test<Type>(w);

        w = (T & w);
    }
}

int main(int argc, char *argv[])
{
    arguments::addBoolOption("parallel", "run in parallel");
    arguments args(argc, argv);

    // // Unit cube

    vertexVector unit;

    for(int i = 0; i < 8; i++)
        unit[i] = unitCube[i];

    // Skewed hex

    vertexVector skewed;

    skewed.lba() = vector(1,   2, 3  );
    skewed.rba() = vector(2,   2, 3  );
    skewed.lta() = vector(1.1, 3, 3.2);
    skewed.rta() = vector(2.1, 3, 3.2);
    skewed.lbf() = vector(1,   2, 4  );
    skewed.rbf() = vector(2,   2, 4  );
    skewed.ltf() = vector(1.1, 3, 4.2);
    skewed.rtf() = vector(2.1, 3, 4.2);

    // General hex

    vertexVector general;

    general.lba() = vector(1,  2,  2.9);
    general.rba() = vector(2.1,2,  2.8);
    general.lta() = vector(1,  3.2,2.7);
    general.rta() = vector(2,  3.3,2.6);
    general.lbf() = vector(0.9,2,  4  );
    general.rbf() = vector(2,  2,  4  );
    general.ltf() = vector(1.2,2.5,4.2);
    general.rtf() = vector(2,  2.9,4  );

    // Stretch tensor

    tensor stretch(2,0,0,0,1.9,0,0,0,0.2);

    // Test rotated hex

    testRotated<hexa>(unit);
    testRotated<hexa>(skewed);
    testRotated<hexa>(0.5*skewed);
    testRotated<hexa>(stretch & general);
    testRotated<hexa>(0.9*general);
    testRotated<hexa>(0.1*general+unit);

    // Test rotated piped

    testRotated<piped>(unit);
    testRotated<piped>(skewed);
    testRotated<piped>(0.5*skewed);
}
