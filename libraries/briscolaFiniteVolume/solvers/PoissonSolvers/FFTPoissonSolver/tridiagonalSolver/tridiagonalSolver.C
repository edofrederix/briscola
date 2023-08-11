#include "tridiagonalSolver.H"
#include "mathematicalConstants.H"
#include "FFTPoissonSolver.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

namespace FFT
{

using Foam::sin;
using Foam::sqr;

void tridiagonalSolver::computeEigenvalues()
{
    const scalar pi(constant::mathematical::pi);

    lambda1_ = scalarList(Nd_[dir1_], Zero);
    lambda2_ = scalarList(Nd_[dir2_], Zero);

    label Sd1 = Sd_[dir1_];
    label Sd2 = Sd_[dir2_];

    label N1 = N_[dir1_];
    label N2 = N_[dir2_];

    // Set eingenvalues for transform in first direction

    switch (BC_[dir1_])
    {
        case 0:
            break;
        case 1:
            for (int i = 0; i < Nd_[dir1_]; i++)
            {
                lambda1_[i] = - 4.0
                * sqr(sin((Sd1+i+1) * pi / (2.0 * N1)));
            }
            break;

        case 2:
            for (int i = 0; i < Nd_[dir1_]; i++)
            {
                lambda1_[i] = - 4.0
                * sqr(sin((Sd1+i) * pi / (2.0 * N1)));
            }
            break;

        case 3:
        case 4:
            for (int i = 0; i < Nd_[dir1_]; i++)
            {
                lambda1_[i] = - 4.0
                * sqr(sin((2.0*(Sd1+i+1)-1.0) * pi / (4.0 * N1)));
            }
            break;

        case 5:
            for (int i = 0; i < Nd_[dir1_]; i++)
            {
                lambda1_[i] = -4.0
                * sqr(sin((Sd1+i) * pi / N1));
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }

    // Set eingenvalues for transform in second direction

    switch (BC_[dir2_])
    {
        case 0:
            break;
        case 1:
            for (int j = 0; j < Nd_[dir2_]; j++)
            {
                lambda2_[j] = - 4.0
                * sqr(sin((Sd2+j+1) * pi / (2.0 * N2)));
            }
            break;

        case 2:
            for (int j = 0; j < Nd_[dir2_]; j++)
            {
                lambda2_[j] = - 4.0
                * sqr(sin((Sd2+j) * pi / (2.0 * N2)));
            }
            break;

        case 3:
        case 4:
            for (int j = 0; j < Nd_[dir2_]; j++)
            {
                lambda2_[j] = - 4.0
                * sqr(sin((2.0*(Sd2+j+1)-1.0) * pi / (4.0 * N2)));
            }
            break;

        case 5:
            for (int j = 0; j < Nd_[dir2_]; j++)
            {
                lambda2_[j] = -4.0
                * sqr(sin((Sd2+j) * pi / N2));
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }
}

void tridiagonalSolver::computeDiagonals()
{
    // Square of cell sizes in each direction
    scalar d1sqr = sqr(cellSizes_[dir1_][0]);
    scalar d2sqr = sqr(cellSizes_[dir2_][0]);
    scalarList d3sqr = sqr(cellSizes_[solver_.FFTPlan().solveDir()]);

    label Nsolve = N_[solver_.FFTPlan().solveDir()];

    label cursor = 0;

    for (int i = 0; i < Nd_[dir1_]; i++)
    {
        for (int j = 0; j < Nd_[dir2_]; j++)
        {
            // Set values of first and last coefficients
            // of main diagonal based on boundary conditions

            switch (BC_[solver_.FFTPlan().solveDir()])
            {
                case 1:
                    D_(cursor) = lambda1_[i]/d1sqr
                        + lambda2_[j]/d2sqr - 3.0/d3sqr[0];

                    D_(cursor + Nsolve - 1) = lambda1_[i]/d1sqr
                        + lambda2_[j]/d2sqr - 3.0/d3sqr[Nsolve - 1];
                    break;

                case 2:
                    D_(cursor) = lambda1_[i]/d1sqr
                        + lambda2_[j]/d2sqr - 1.0/d3sqr[0];

                    D_(cursor + Nsolve - 1) = lambda1_[i]/d1sqr
                        + lambda2_[j]/d2sqr - 1.0/d3sqr[Nsolve - 1];

                    // If coefficient matrix is singular
                    if (lambda1_[i] == 0 && lambda2_[j] == 0)
                    {
                        D_(cursor) -= 1e-10;
                    }
                    break;

                case 3:
                    D_(cursor) = lambda1_[i]/d1sqr
                        + lambda2_[j]/d2sqr - 3.0/d3sqr[0];

                    D_(cursor + Nsolve - 1) = lambda1_[i]/d1sqr
                        + lambda2_[j]/d2sqr - 1.0/d3sqr[Nsolve - 1];
                    break;

                case 4:
                    D_(cursor) = lambda1_[i]/d1sqr
                        + lambda2_[j]/d2sqr - 1.0/d3sqr[0];

                    D_(cursor + Nsolve - 1) = lambda1_[i]/d1sqr
                        + lambda2_[j]/d2sqr - 3.0/d3sqr[Nsolve - 1];
                    break;

                case 5:
                    D_(cursor) = lambda1_[i]/d1sqr
                        + lambda2_[j]/d2sqr - 2.0/d3sqr[0];
                    D_(cursor + Nsolve - 1) = lambda1_[i]/d1sqr
                        + lambda2_[j]/d2sqr - 2.0/d3sqr[Nsolve - 1];

                    // If coefficient matrix is singular
                    if (lambda1_[i] == 0 && lambda2_[j] == 0)
                    {
                        D_(cursor) -= 1e-10;
                    }
                    break;

                default:
                    FatalError
                        << "Incorrect solve direction." << endl
                        << abort(FatalError);
                    break;
            }

            // Set remaining main diagonal values

            for (int k = 1; k < Nsolve - 1; k++)
            {
                D_(cursor + k) = lambda1_[i]/d1sqr
                    + lambda2_[j]/d2sqr - 2.0/d3sqr[k];
            }

            cursor += Nsolve;
        }
    }
}

tridiagonalSolver::tridiagonalSolver
(
    FFTPoissonSolver& solver,
    const labelVector& BC
)
:
    solver_(solver),
    N_(solver.fvMsh().msh().cast<rectilinearMesh>().N()),
    Nd_(solver.decomp().Nd(solver_.FFTPlan().solveDir())[Pstream::myProcNo()]),
    Sd_(solver.decomp().Sd(solver_.FFTPlan().solveDir())[Pstream::myProcNo()]),
    BC_(BC),
    cellSizes_(solver.fvMsh().msh().cast<rectilinearMesh>().globalCellSizes()),
    D_(Nd_, Zero),
    DU_(1.0 / sqr(cellSizes_[solver_.FFTPlan().solveDir()])),
    DL_(1.0 / sqr(cellSizes_[solver_.FFTPlan().solveDir()]))
{
    switch (solver_.FFTPlan().solveDir())
    {
        case 0:
            dir1_ = 2;
            dir2_ = 1;
            break;

        case 1:
            dir1_ = 0;
            dir2_ = 2;
            break;

        case 2:
            dir1_ = 0;
            dir2_ = 1;
            break;

        default:
            FatalError
                << "Incorrect solve direction." << endl
                << abort(FatalError);
            break;
    }

    computeEigenvalues();
    computeDiagonals();
}

tridiagonalSolver::~tridiagonalSolver()
{}

void tridiagonalSolver::solve
(
    scalarBlock& p,
    const scalarBlock& f,
    const bool ddt,
    const label start
)
{
    const scalar deltaT = solver_.fvMsh().time().deltaTValue();

    labelVector N(p.shape());

    label Nsolve = N_[solver_.FFTPlan().solveDir()];

    label cursor = 0;

    for (int i = 0; i < Nd_[dir1_]; i++)
    {
        for (int j = 0; j < Nd_[dir2_]; j++)
        {
            // Copy D and f so as not to overwrite them
            scalarList Dh(Nsolve, Zero);
            scalarList fh(Nsolve, Zero);

            for (int k = 0; k < Nsolve; k++)
            {
                Dh[k] = D_(cursor + k) - (ddt ? 1.0/deltaT : 0.0);
                fh[k] = f(cursor + k);
            }

            // Forward substitution
            for (int k = 1 + start; k < Nsolve; k++)
            {
                scalar m = DL_[k] / Dh[k-1];
                Dh[k] -= m * DU_[k-1];
                fh[k] -= m * fh[k-1];
            }

            // Backward substitution
            p(cursor + Nsolve-1) = fh[Nsolve-1] / Dh[Nsolve-1];

            for (int k = Nsolve-2; k >= start; k--)
            {
                p(cursor + k) = (fh[k] - DU_[k] * p(cursor + k+1)) / Dh[k];
            }

            // The cursor moves along the two unsolved directions
            cursor += Nsolve;
        }
    }
}

void tridiagonalSolver::solveCyclic
(
    scalarBlock& p,
    const scalarBlock& f,
    const bool ddt
)
{
    const scalar deltaT = solver_.fvMsh().time().deltaTValue();

    // Solve first auxiliary system
    scalarBlock u(Nd_, Zero);

    solve(u, f, ddt, 1);

    // Solve second auxiliary system
    scalarBlock v(Nd_, Zero);
    scalarBlock vf(Nd_, Zero);

    const label Nsolve = N_[solver_.FFTPlan().solveDir()];

    label cursor = 0;

    for (int i = 0; i < Nd_[dir1_]; i++)
    {
        for (int j = 0; j < Nd_[dir2_]; j++)
        {
            vf(cursor + 1) = -DL_[1];
            vf(cursor + Nsolve-1) = -DU_[Nsolve-1];

            cursor += Nsolve;
        }
    }

    solve(v, vf, 1);

    cursor = 0;

    // Reconstruct solution
    for (int i = 0; i < Nd_[dir1_]; i++)
    {
        for (int j = 0; j < Nd_[dir2_]; j++)
        {
            p(cursor) =
                (
                      f(cursor)
                    - DL_[0] * u(cursor + Nsolve-1)
                    - DU_[0] * u(cursor + 1)
                )
                /
                (
                      D_(cursor) - (ddt ? 1.0/deltaT : 0.0)
                    + DL_[0] * v(cursor + Nsolve-1)
                    + DU_[0] * v(cursor + 1)
                );

            for (int k = 1; k < Nsolve; k++)
            {
                p(cursor + k) = u(cursor + k) + p(cursor) * v(cursor + k);
            }

            cursor += Nsolve;
        }
    }
}

void tridiagonalSolver::solve(scalarBlock& xyzPencil, const bool ddt)
{
    // Copy RHS of tridiagonal system
    scalarBlock f(xyzPencil);

    if (BC_[solver_.FFTPlan().solveDir()] != 5)
    {
        solve(xyzPencil, f, ddt);
    }
    else
    {
        solveCyclic(xyzPencil, f, ddt);
    }
}

}

}

}

}