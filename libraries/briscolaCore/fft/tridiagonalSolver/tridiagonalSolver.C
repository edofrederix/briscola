#include "tridiagonalSolver.H"

namespace Foam
{

namespace briscola
{

// Constructors

tridiagonalSolver::tridiagonalSolver(label n)
:
    N_(n)
{}

// Destructor

tridiagonalSolver::~tridiagonalSolver()
{}

// Public

void tridiagonalSolver::tridiagonalSolve(label n, double DU[], double D[], double DL[], double f[], double p[])
{
    List<scalar> Dh(n);
    List<scalar> fh(n);

    // Copy D and f arrays so as not to overwrite them
    for(int i = 0; i < n; i++)
    {
        Dh[i] = D[i];
        fh[i] = f[i];
    }

    // Forward substitution
    for (int i = 1; i < n; i++)
    {
        scalar m = DL[i] / Dh[i-1];
        Dh[i] -= m * DU[i-1];
        fh[i] -= m * fh[i-1];
    }

    // Backward substitution
    p[n-1] = fh[n-1] / Dh[n-1];
    for (int i = n-2; i >= 0; i--)
    {
        p[i] = (fh[i] - DU[i] * p[i+1]) / Dh[i];
    }
}

void tridiagonalSolver::cyclicTridiagonalSolve(label n, double DU[], double D[], double DL[], double f[], double p[])
{
    List<scalar> u(n, Zero);
    tridiagonalSolve((n-1), &DU[1], &D[1], &DL[1], &f[1], &u[1]);

    List<scalar> v(n, Zero);
    List<scalar> vf(n, Zero);

    vf[1] = -DL[1];
    vf[n-1] = -DU[n-1];

    tridiagonalSolve((n-1), &DU[1], &D[1], &DL[1], &vf[1], &v[1]);

    p[0] = (f[0] - DL[0] * u[n-1] - DU[0] * u[1]) / (D[0] + DL[0] * v[n-1] + DU[0] * v[1]);

    for (int k = 1; k < n; k++)
    {
        p[k] = u[k] + p[0]*v[k];
    }
}

}

}
