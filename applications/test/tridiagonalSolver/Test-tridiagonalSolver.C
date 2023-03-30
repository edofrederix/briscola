#include "arguments.H"
#include "Time.H"
#include "vector.H"
#include "block.H"
#include "fv.H"
#include "tridiagonalSolver.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    using std::rand;
    using std::srand;

    arguments args(argc, argv);

    // Create case and construct tridiagonal solver

    #include "createBriscolaTime.H"

    labelVector N(64,64,64);

    labelVector Sz(0,0,0);
    labelVector BC(1,1,1);

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

    tridiagonalSolver tds(fvMsh, N, Sz, BC);

    // Input data containing random values

    int seed = 123;
    srand(seed);

    scalarList DU(N.z());
    scalarList DL(N.z());

    scalarBlock D(N);
    scalarBlock f(N);

    for (int i = 0; i < N.x(); i++)
    {
        for (int j = 0; j < N.y(); j++)
        {
            for (int k = 0; k < N.z(); k++)
            {
                D(i,j,k) = static_cast<double>(rand()) / RAND_MAX;
                f(i,j,k) = static_cast<double>(rand()) / RAND_MAX;
                DU[k] = static_cast<double>(rand()) / RAND_MAX;
                DL[k] = static_cast<double>(rand()) / RAND_MAX;
            }
        }
    }

    scalarBlock p(N);

    // Test tridiagonal solver

    tds.solve
    (
        DU,
        D,
        DL,
        f,
        p
    );

    for (int i = 0; i < N.x(); i++)
    {
        for (int j = 0; j < N.y(); j++)
        {
            for (int k = 0; k < N.z(); k++)
            {
                scalar residual =
                  ((k == 0)       ? 0 : DL[k] * p(i,j,k-1))
                + ((k == N.z()-1) ? 0 : DU[k] * p(i,j,k+1))
                + (D(i,j,k) * p(i,j,k) - f(i,j,k));

                if(mag(residual) > 1e-5)
                {
                    FatalError
                        << "Test failed. Residual =  " << residual
                        << " at index " << labelVector(i,j,k) << endl
                        << abort(FatalError);
                }
            }
        }
    }

    Info << "--------------------------------------------" << nl;
    Info << "Tridiagonal system solution check successful" << nl;
    Info << "--------------------------------------------" << endl;

    // Test cyclic tridiagonal solver

    tds.solveCyclic
    (
        DU,
        D,
        DL,
        f,
        p
    );

    Info << nl;

    for (int i = 0; i < N.x(); i++)
    {
        for (int j = 0; j < N.y(); j++)
        {
            for (int k = 0; k < N.z(); k++)
            {
                scalar residual =
                  ((k == 0)       ? DL[k] * p(i,j,N.z()-1) : DL[k] * p(i,j,k-1))
                + ((k == N.z()-1) ? DU[N.z()-1] * p(i,j,0) : DU[k] * p(i,j,k+1))
                + (D(i,j,k) * p(i,j,k) - f(i,j,k));

                if(mag(residual) > 1e-5)
                {
                    FatalError
                        << "Test failed. Residual =  " << residual
                        << " at index " << labelVector(i,j,k) << endl
                        << abort(FatalError);
                }
            }
        }
    }

    Info << "---------------------------------------------------" << nl;
    Info << "Cyclic tridiagonal system solution check successful" << nl;
    Info << "---------------------------------------------------" << endl;
}
