#include "arguments.H"
#include "Time.H"
#include "fvMesh.H"
#include "vof.H"

#include "truncatedPiped.H"
#include "truncatedHex.H"
#include "SortableList.H"
#include "constants.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

void testVolume(const vertexVector& v, const vector n)
{
    SortableList<scalar> V(8);
    scalarList C(8);

    for (int i = 0; i < 8; i++)
    {
        C[i] = - (n & v[i]);
        V[i] = truncatedPiped(v,n,C[i]).volume();
    }

    V.sort();

    // Test for monotonicity of the solution in the C brackets

    scalar Vi = -1;
    const label N = 20;

    for (int i = 1; i < 8; i++)
    {
        scalar CMin = C[V.indices()[i-1]];
        scalar CMax = C[V.indices()[i]];

        for (int j = 0; j <= N; j++)
        {
            const scalar C = CMin + (CMax-CMin)*j/N;

            truncatedPiped thex(v,n,C);

            const scalar VNew = thex.volume();

            if (VNew - Vi < -1e-12)
                FatalErrorInFunction
                    << "Test 1 failed" << endl << abort(FatalError);

            Vi = VNew;
        }
    }
}

void testRotatedVolumes(const vertexVector& v, const vector n)
{
    const scalar theta = 5.0/360*2*Foam::constant::mathematical::pi;

    // Rotate around x

    tensor T
    (
        vector(1,0,0),
        vector(0,Foam::cos(theta),-Foam::sin(theta)),
        vector(0,Foam::sin(theta),Foam::cos(theta))
    );

    vertexVector hex = v;

    for (int i = 0; i < 72; i++)
    {
        testVolume(hex, n);

        hex = (T & hex);
    }

    // Rotate around y

    T = tensor
    (
        vector(Foam::cos(theta),0,Foam::sin(theta)),
        vector(0,1,0),
        vector(-Foam::sin(theta),0,Foam::cos(theta))
    );

    hex = v;

    for (int i = 0; i < 72; i++)
    {
        testVolume(hex, n);

        hex = (T & hex);
    }

    // Rotate around z

    T = tensor
    (
        vector(Foam::cos(theta),-Foam::sin(theta),0),
        vector(Foam::sin(theta),Foam::cos(theta),0),
        vector(0,0,1)
    );

    hex = v;

    for (int i = 0; i < 72; i++)
    {
        testVolume(hex, n);

        hex = (T & hex);
    }
}

void testLVE
(
    const vof& vf,
    const vertexVector& v,
    const vector n,
    const word method
)
{
    if (method != "p" && method != "a")
        FatalErrorInFunction
            << "Invalid LVE method. Should be p or a."
            << endl << abort(FatalError);

    SortableList<scalar> C(8);
    scalarList V(8);

    for (int i = 0; i < 8; i++)
    {
        C[i] = - (n & v[i]);
        V[i] = truncatedHex(v,n,C[i]).volume();
    }

    C.sort();

    const scalar Vc = Foam::max(V);

    // Test inverse solution exactly at vertices

    for (int i = 1; i < 8; i++)
    {
        const scalar fi = V[C.indices()[i]]/Vc;

        const scalar Ci =
            method == "p" ? vf.lve().pLVE(v,n,fi)
                          : vf.lve().aLVE(v,n,fi);

        const scalar fj = truncatedHex(v,n,Ci).volume()/Vc;

        if (Foam::mag(fi - fj) > 1e-8)
            FatalErrorInFunction
                << fi << " " << fj << " Test 2a failed" << endl << abort(FatalError);
    }

    // Test inverse solution over C range

    const scalar CMin = C[0];
    const scalar CMax = C[7];

    int N = 20;
    for (int i = 0; i <= N; i++)
    {
        const scalar C1 = CMin + (CMax-CMin)*i/N;
        const scalar f1 = truncatedHex(v,n,C1).volume()/Vc;

        const scalar C2 =
            method == "p" ? vf.lve().pLVE(v,n,f1)
                          : vf.lve().aLVE(v,n,f1);

        const scalar f2 = truncatedHex(v,n,C2).volume()/Vc;

        if (Foam::mag(f1-f2) > 1e-8)
            FatalErrorInFunction
                << "Test 2b failed" << endl << abort(FatalError);
    }

    // Check if sign is correct

    const scalar Cf =
        method == "p" ? vf.lve().pLVE(v,n,1)
                      : vf.lve().aLVE(v,n,1);

    if (Foam::mag(truncatedHex(v,n,Cf).volume()/Vc - 1.0) > 1e-8)
        FatalErrorInFunction
                << "Test 2c failed" << endl << abort(FatalError);

    const scalar Ce =
        method == "p" ? vf.lve().pLVE(v,n,0)
                      : vf.lve().aLVE(v,n,0);

    if (Foam::mag(truncatedHex(v,n,Ce).volume()/Vc) > 1e-8)
        FatalErrorInFunction
                << "Test 2d failed" << endl << abort(FatalError);
}

void testRotatedLVE
(
    const vof& vf,
    const vertexVector& v,
    const vector n,
    const word method
)
{
    const scalar theta = 5.0/360*2*Foam::constant::mathematical::pi;

    // Rotate around x

    tensor T
    (
        vector(1,0,0),
        vector(0,Foam::cos(theta),-Foam::sin(theta)),
        vector(0,Foam::sin(theta),Foam::cos(theta))
    );

    vertexVector hex = v;

    for (int i = 0; i < 72; i++)
    {
        testLVE(vf, hex, n, method);

        hex = (T & hex);
    }

    // Rotate around y

    T = tensor
    (
        vector(Foam::cos(theta),0,Foam::sin(theta)),
        vector(0,1,0),
        vector(-Foam::sin(theta),0,Foam::cos(theta))
    );

    hex = v;

    for (int i = 0; i < 72; i++)
    {
        testLVE(vf, hex, n, method);

        hex = (T & hex);
    }

    // Rotate around z

    T = tensor
    (
        vector(Foam::cos(theta),-Foam::sin(theta),0),
        vector(Foam::sin(theta),Foam::cos(theta),0),
        vector(0,0,1)
    );

    hex = v;

    for (int i = 0; i < 72; i++)
    {
        testLVE(vf, hex, n, method);

        hex = (T & hex);
    }
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createBriscolaVof.H"

    // Unit cube

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

    // Stretch tensor

    tensor stretch(2,0,0,0,1.9,0,0,0,0.2);

    forAll(ns, ni)
    {
        const vector n = ns[ni];

        // Volume computation

        testRotatedVolumes(unit, n);
        testRotatedVolumes(skewed, n);
        testRotatedVolumes(0.5*skewed, n);
        testRotatedVolumes(stretch & general, n);
        testRotatedVolumes(0.9*general, n);
        testRotatedVolumes(0.1*general+unit, n);

        // Local volume enforcement

        // Analytical algorithm of Scardovelli & Zaleski (works only on
        // parallelepipeds)

        /*

        testRotatedLVE(vf, unit, n, "p");
        testRotatedLVE(vf, 0.5*unit, n, "p");
        testRotatedLVE(vf, tensor(1.1,0,0,0,1,0,0,0,1) & unit, n, "p");
        testRotatedLVE(vf, stretch & unit, n, "p");
        testRotatedLVE(vf, unit+unit, n, "p");
        testRotatedLVE(vf, skewed, n, "p");
        testRotatedLVE(vf, skewed + unit, n, "p");
        testRotatedLVE(vf, stretch & (skewed + unit), n, "p");

        // General algorithm for arbitrary hexahedrons

        testRotatedLVE(vf, 0.6*unit, n, "a");
        testRotatedLVE(vf, 0.9*skewed, n, "a");
        testRotatedLVE(vf, 0.45*general, n, "a");
        testRotatedLVE(vf, 0.45*(general+unit), n, "a");
        testRotatedLVE(vf, 0.45*(general-unit), n, "a");

        // Test cells

        label N = 100;

        const colocatedVertexVectorDirection& v =
            fvMsh.template metrics<colocated>().vertexCenters()[0][0];

        const colocatedScalarDirection& V =
            fvMsh.template metrics<colocated>().cellVolumes()[0][0];

        for (int i = 0; i <= N; i++)
        {
            colocatedScalarDirection& a = vf.alpha()[0][0];

            a = scalar(i)/N;

            forAllCells(a, i, j, k)
            {
                const scalar C = vf.lve()(i,j,k,n);

                const scalar f = truncatedHex(v(i,j,k),n,C).volume()/V(i,j,k);

                if (Foam::mag(f-a(i,j,k)) > 1e-8)
                    FatalErrorInFunction
                        << "Test 3 failed" << endl << abort(FatalError);
            }
        }

        */
    }
}