#include "arguments.H"
#include "vector.H"
#include "block.H"
#include "tridiagonalSolver.H"

using namespace Foam;
using namespace briscola;

int main(int argc, char *argv[])
{
    using std::rand;
    using std::srand;

    arguments::validArgs.append("size");

    arguments args(argc, argv);

    const label N(args.argRead<label>(1));


    // Input data containing random values

    int seed = 123;
    srand(seed);

    List<scalar> DU(N);
    List<scalar> D (N);
    List<scalar> DL(N);
    List<scalar> f (N);

    List<scalar> p (N);


    for (int i = 0; i < N; i++)
    {
        DU[i] = static_cast<double>(rand()) / RAND_MAX;
        D [i] = static_cast<double>(rand()) / RAND_MAX;
        DL[i] = static_cast<double>(rand()) / RAND_MAX;
        f [i] = static_cast<double>(rand()) / RAND_MAX;
    }

    autoPtr<tridiagonalSolver> tds;
    tds = new tridiagonalSolver(N);

    tds->tridiagonalSolve
    (
        N,
        reinterpret_cast<double*>(DU.begin()),
        reinterpret_cast<double*>(D .begin()),
        reinterpret_cast<double*>(DL.begin()),
        reinterpret_cast<double*>(f .begin()),
        reinterpret_cast<double*>(p .begin())
    );

    for (int i = 0; i < N; i++)
    {
        scalar residual =
          ((i == 0)   ?          0 : DL[i] * p[i-1])
        + ((i == N-1) ?          0 : DU[i] * p[i+1])
        + (D[i] * p[i] - f[i]);

        if(mag(residual) > 1e-10)
        {
            FatalError
                << "Test failed. Residual =  " << residual
                << " at index " << i << endl
                << abort(FatalError);
        }
    }

    Info << "--------------------------------------------" << nl;
    Info << "Tridiagonal system solution check successful" << nl;
    Info << "--------------------------------------------" << endl;


    tds->cyclicTridiagonalSolve
    (
        N,
        reinterpret_cast<double*>(DU.begin()),
        reinterpret_cast<double*>(D .begin()),
        reinterpret_cast<double*>(DL.begin()),
        reinterpret_cast<double*>(f .begin()),
        reinterpret_cast<double*>(p .begin())
    );

    Info << nl;

    for (int i = 0; i < N; i++)
    {
        scalar residual =
          ((i == 0)   ?          DL[i] * p[N-1] : DL[i] * p[i-1])
        + ((i == N-1) ?          DU[N-1] * p[0] : DU[i] * p[i+1])
        + (D[i] * p[i] - f[i]);

        if(mag(residual) > 1e-10)
        {
            FatalError
                << "Test failed. Residual =  " << residual
                << " at index " << i << endl
                << abort(FatalError);
        }
    }

    Info << "---------------------------------------------------" << nl;
    Info << "Cyclic tridiagonal system solution check successful" << nl;
    Info << "---------------------------------------------------" << endl;
}
