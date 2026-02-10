#include "arguments.H"
#include "Time.H"
#include "fvMesh.H"
#include "twoPhaseModel.H"
#include "geometricVof.H"
#include "geometricObjects.H"
#include "SortableList.H"
#include "constants.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

const label Ntheta = 8;
const label NC = 8;

void testVolumeHex(const vertexVector& v, const vector n)
{
    SortableList<scalar> V(8);
    scalarList C(8);

    for (int i = 0; i < 8; i++)
    {
        C[i] = - (n & v[i]);
        V[i] = hexa(v).truncationVolume(n,C[i]);
    }

    V.sort();

    // Test for monotonicity of the solution in the C brackets

    scalar Vi = -1;

    for (int i = 1; i < 8; i++)
    {
        scalar CMin = C[V.indices()[i-1]];
        scalar CMax = C[V.indices()[i]];

        for (int j = 0; j <= NC; j++)
        {
            const scalar C = CMin + (CMax-CMin)*j/scalar(NC);

            const scalar VNew = hexa(v).truncationVolume(n,C);

            if (VNew - Vi < -1e-12)
                FatalErrorInFunction
                    << "Test 1 failed" << endl << abort(FatalError);

            Vi = VNew;
        }
    }
}

void testVolumePiped(const vertexVector& v, const vector n)
{
    SortableList<scalar> V(8);
    scalarList C(8);

    for (int i = 0; i < 8; i++)
    {
        C[i] = - (n & v[i]);
        V[i] = piped(v).truncationVolume(n,C[i]);
    }

    V.sort();

    // Test for monotonicity of the solution in the C brackets

    scalar Vi = -1;

    for (int i = 1; i < 8; i++)
    {
        scalar CMin = C[V.indices()[i-1]];
        scalar CMax = C[V.indices()[i]];

        for (int j = 0; j <= NC; j++)
        {
            const scalar C = CMin + (CMax-CMin)*j/scalar(NC);

            const scalar VNew = piped(v).truncationVolume(n,C);

            if (VNew - Vi < -1e-12)
                FatalErrorInFunction
                    << "Test 1 failed" << endl << abort(FatalError);

            Vi = VNew;
        }
    }
}

void testRotatedVolumesHex(const vertexVector& v, const vector n)
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

    vertexVector hex = v;

    for (int i = 0; i < Ntheta; i++)
    {
        testVolumeHex(hex, n);

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

    for (int i = 0; i < Ntheta; i++)
    {
        testVolumeHex(hex, n);

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

    for (int i = 0; i < Ntheta; i++)
    {
        testVolumeHex(hex, n);

        hex = (T & hex);
    }
}

void testRotatedVolumesPiped(const vertexVector& v, const vector n)
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

    vertexVector hex = v;

    for (int i = 0; i < Ntheta; i++)
    {
        testVolumePiped(hex, n);

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

    for (int i = 0; i < Ntheta; i++)
    {
        testVolumePiped(hex, n);

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

    for (int i = 0; i < Ntheta; i++)
    {
        testVolumePiped(hex, n);

        hex = (T & hex);
    }
}

void testLVE
(
    const geometricVof& gvf,
    const vertexVector& v,
    const vector n,
    const word method
)
{

    if (method != "p" && method != "c")
        FatalErrorInFunction
            << "Invalid LVE method. Should be p or c."
            << endl << abort(FatalError);

    SortableList<scalar> C(8);
    scalarList V(8);

    for (int i = 0; i < 8; i++)
    {
        C[i] = - (n & v[i]);
        V[i] = hexa(v).truncationVolume(n,C[i]);
    }

    C.sort();

    const scalar Vc = Foam::max(V);

    // Test inverse solution exactly at vertices

    for (int i = 1; i < 8; i++)
    {
        const scalar fi = V[C.indices()[i]]/Vc;
        scalar Ci;

        if (method == "p")
        {
            Ci = gvf.lve().pLVE(v,n,fi);
        }
        else
        {
            Ci = gvf.lve().cLVE(v,n,fi);
        }

        const scalar fj = hexa(v).truncationVolume(n,Ci)/Vc;

        if (Foam::mag(fi - fj) > 1e-7)
            FatalErrorInFunction
                << fi << " " << fj << " Test 2a failed" << endl << abort(FatalError);
    }

    // Test inverse solution over C range

    const scalar CMin = C[0];
    const scalar CMax = C[7];

    int N = 20;
    for (int i = 0 ; i <= N; i++)
    {
        const scalar C1 = CMin + 0.3 * (CMax-CMin)*i/N;
        const scalar f1 = hexa(v).truncationVolume(n,C1)/Vc;
        scalar C2;

        if (method == "p")
        {
            C2 = gvf.lve().pLVE(v,n,f1);
        }
        else
        {
            C2 = gvf.lve().cLVE(v,n,f1);
        }

        const scalar f2 = hexa(v).truncationVolume(n,C2)/Vc;

        if (Foam::mag(f1-f2) > 1e-7)
            FatalErrorInFunction
                << "Test 2b failed" << endl << abort(FatalError);
    }

    // Check if sign is correct

    scalar Cf, Ce;
    if (method == "p")
    {
        Cf = gvf.lve().pLVE(v,n,1);
    }
    else
    {
        Cf = gvf.lve().cLVE(v,n,1);
    }

    if (Foam::mag(hexa(v).truncationVolume(n,Cf)/Vc - 1.0) > 1e-7)
        FatalErrorInFunction
                << "Test 2c failed" << endl << abort(FatalError);

    if (method == "p")
    {
        Ce = gvf.lve().pLVE(v,n,0);
    }
    else
    {
        Ce = gvf.lve().cLVE(v,n,0);
    }

    if (Foam::mag(hexa(v).truncationVolume(n,Ce)/Vc) > 1e-7)
        FatalErrorInFunction
                << "Test 2d failed" << endl << abort(FatalError);

}

void testRotatedLVE
(
    const geometricVof& gvf,
    const vertexVector& v,
    const vector n,
    const word method
)
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

    vertexVector hex = v;

    for (int i = 0; i < Ntheta; i++)
    {
        testLVE(gvf, hex, n, method);

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

    for (int i = 0; i < Ntheta; i++)
    {
        testLVE(gvf, hex, n, method);

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

    for (int i = 0; i < Ntheta; i++)
    {
        testLVE(gvf, hex, n, method);

        hex = (T & hex);
    }
}

int main(int argc, char *argv[])
{
    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"
    #include "createBriscolaMesh.H"
    #include "createBriscolaVof.H"

    // Vof object must be castable to geometricVof

    geometricVof& gvf = vf.cast<geometricVof>();

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

        // Volume computation (TruncatedHex)

        testRotatedVolumesHex(unit, n);
        testRotatedVolumesHex(skewed, n);
        testRotatedVolumesHex(0.5*skewed, n);
        testRotatedVolumesHex(stretch & general, n);
        testRotatedVolumesHex(0.9*general, n);
        testRotatedVolumesHex(0.1*general+unit, n);

        // Volume computation (TruncatedPiped)

        testRotatedVolumesPiped(unit, n);
        testRotatedVolumesPiped(skewed, n);
        testRotatedVolumesPiped(0.5*skewed, n);

        // Local volume enforcement

        // Analytical algorithm of Scardovelli & Zaleski (works only on
        // parallelepipeds)

        testRotatedLVE(gvf, unit, n, "p");
        testRotatedLVE(gvf, 0.5*unit, n, "p");
        testRotatedLVE(gvf, tensor(1.1,0,0,0,1,0,0,0,1) & unit, n, "p");
        testRotatedLVE(gvf, stretch & unit, n, "p");
        testRotatedLVE(gvf, unit+unit, n, "p");
        testRotatedLVE(gvf, skewed, n, "p");
        testRotatedLVE(gvf, skewed + unit, n, "p");
        testRotatedLVE(gvf, stretch & (skewed + unit), n, "p");

        // General algorithm for arbitrary hexahedrons

        testRotatedLVE(gvf, 0.6*unit, n, "c");
        testRotatedLVE(gvf, 0.9*skewed, n, "c");
        testRotatedLVE(gvf, 0.45*general, n, "c");
        testRotatedLVE(gvf, 0.45*(general+unit), n, "c");

        // Test cells

        const colocatedVertexVectorField& v =
            fvMsh.template metrics<colocated>().vertexCenters();

        const colocatedScalarField& cv =
            fvMsh.template metrics<colocated>().cellVolumes();

        for (int i = 0; i <= NC; i++)
        {
            colocatedScalarField& a = gvf.alpha();

            a = scalar(i)/NC;

            forAllCells(a, i, j, k)
            {
                const scalar C = gvf.lve()(a(i,j,k),v(i,j,k),n);

                const scalar f = hexa(v(i,j,k)).truncationVolume(n,C)/cv(i,j,k);

                if (Foam::mag(f-a(i,j,k)) > 1e-8)
                    FatalErrorInFunction
                        << "Test 3 failed" << endl << abort(FatalError);
            }
        }
    }
}
