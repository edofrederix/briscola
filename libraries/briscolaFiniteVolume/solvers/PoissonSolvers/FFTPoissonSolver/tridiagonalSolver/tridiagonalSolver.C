#include "tridiagonalSolver.H"
#include "mathematicalConstants.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

void tridiagonalSolver::computeEigenvalues()
{
    using Foam::sin;

    const scalar pi(constant::mathematical::pi);

    label dir1 = 0;
    label dir2 = 0;

    switch (solveDir_)
    {
        case 0:
            dir1 = 2;
            dir2 = 1;
            break;

        case 1:
            dir1 = 0;
            dir2 = 2;
            break;

        case 2:
            dir1 = 0;
            dir2 = 1;
            break;
    }

    label N1 = N_[dir1];
    label N2 = N_[dir2];

    label Nd1 = Nd_[dir1];
    label Nd2 = Nd_[dir2];

    label Sd1 = Sd_[dir1];
    label Sd2 = Sd_[dir2];

    lambda1_ = scalarList(Nd1, Zero);
    lambda2_ = scalarList(Nd2, Zero);

    switch (BC_[dir1])
    {
        case 1:
            for (int i = 0; i < Nd1; i++)
            {
                lambda1_[i] = - 4.0
                * sqr(sin((Sd1+i+1) * pi / (2.0 * N1)));
            }
            break;

        case 2:
            for (int i = 0; i < Nd1; i++)
            {
                lambda1_[i] = - 4.0
                * sqr(sin((Sd1+i) * pi / (2.0 * N1)));
            }
            break;

        case 3:
        case 4:
            for (int i = 0; i < Nd1; i++)
            {
                lambda1_[i] = - 4.0
                * sqr(sin((2.0*(Sd1+i+1)-1.0) * pi / (4.0 * N1)));
            }
            break;

        case 5:
            for (int i = 0; i < Nd1; i++)
            {
                lambda1_[i] = -4.0
                * sqr(sin((Sd1+i) * pi / N1));
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }

    switch (BC_[dir2])
    {
        case 1:
            for (int j = 0; j < Nd2; j++)
            {
                lambda2_[j] = - 4.0
                * sqr(sin((Sd2+j+1) * pi / (2.0 * N2)));
            }
            break;

        case 2:
            for (int j = 0; j < Nd2; j++)
            {
                lambda2_[j] = - 4.0
                * sqr(sin((Sd2+j) * pi / (2.0 * N2)));
            }
            break;

        case 3:
        case 4:
            for (int j = 0; j < Nd2; j++)
            {
                lambda2_[j] = - 4.0
                * sqr(sin((2.0*(Sd2+j+1)-1.0) * pi / (4.0 * N2)));
            }
            break;

        case 5:
            for (int j = 0; j < Nd2; j++)
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
    scalar d1sqr = 0;
    scalar d2sqr = 0;
    scalarList d3sqr = sqr(cellSizes_[solveDir_]);

    switch (solveDir_)
    {
        case 0:
            d1sqr = sqr(cellSizes_[2][0]);
            d2sqr = sqr(cellSizes_[1][0]);
            break;

        case 1:
            d1sqr = sqr(cellSizes_[0][0]);
            d2sqr = sqr(cellSizes_[2][0]);
            break;

        case 2:
            d1sqr = sqr(cellSizes_[0][0]);
            d2sqr = sqr(cellSizes_[1][0]);
            break;
    }

    // Dimension of solved direction
    label Nsolve = 0;

    // Dimensions of unsolved directions
    label N1 = 0;
    label N2 = 0;

    switch (solveDir_)
    {
        case 0:
            Nsolve = Nd_.x();
            N1 = Nd_.z();
            N2 = Nd_.y();
            break;

        case 1:
            Nsolve = Nd_.y();
            N1 = Nd_.x();
            N2 = Nd_.z();
            break;

        case 2:
            Nsolve = Nd_.z();
            N1 = Nd_.x();
            N2 = Nd_.y();
            break;
    }

    label cursor = 0;

    for (int i = 0; i < N1; i++)
    {
        for (int j = 0; j < N2; j++)
        {
            switch (BC_[solveDir_])
            {
                case 1:
                    D_(i,j,0) = lambda1_[i]/d1sqr + lambda2_[j]/d2sqr - 3.0/d3sqr[0];
                    D_(i,j,N_.z()-1) = lambda1_[i]/d1sqr
                                     + lambda2_[j]/d2sqr - 3.0/d3sqr[N_.z()-1];
                    break;

                case 2:
                    D_(i,j,0) = lambda1_[i]/d1sqr + lambda2_[j]/d2sqr - 1.0/d3sqr[0];
                    D_(i,j,N_.z()-1) = lambda1_[i]/d1sqr
                                     + lambda2_[j]/d2sqr - 1.0/d3sqr[N_.z()-1];
                    // If coefficient matrix is singular
                    if (lambda1_[i] == 0 && lambda2_[j] == 0)
                    {
                        D_(i,j,0) -= 1e-10;
                    }
                    break;

                case 3:
                    D_(i,j,0) = lambda1_[i]/d1sqr + lambda2_[j]/d2sqr - 3.0/d3sqr[0];
                    D_(i,j,N_.z()-1) = lambda1_[i]/d1sqr
                                     + lambda2_[j]/d2sqr - 1.0/d3sqr[N_.z()-1];
                    break;

                case 4:
                    D_(i,j,0) = lambda1_[i]/d1sqr + lambda2_[j]/d2sqr - 1.0/d3sqr[0];
                    D_(i,j,N_.z()-1) = lambda1_[i]/d1sqr
                                     + lambda2_[j]/d2sqr - 3.0/d3sqr[N_.z()-1];
                    break;

                case 5:
                    D_(i,j,0) = lambda1_[i]/d1sqr + lambda2_[j]/d2sqr - 2.0/d3sqr[0];
                    D_(i,j,N_.z()-1) = lambda1_[i]/d1sqr
                                     + lambda2_[j]/d2sqr - 2.0/d3sqr[N_.z()-1];
                    // If coefficient matrix is singular
                    if (lambda1_[i] == 0 && lambda2_[j] == 0)
                    {
                        D_(i,j,0) -= 1e-10;
                    }
                    break;
            }

            for (int k = 1; k < Nsolve - 1; k++)
            {
                D_(cursor++) = lambda1_[i]/d1sqr + lambda2_[j]/d2sqr - 2.0/d3sqr[k];
            }
        }
    }
}

tridiagonalSolver::tridiagonalSolver
(
    const fvMesh& fvMsh,
    label solveDir,
    labelVector Nd,
    labelVector Sd,
    labelVector BC
)
:
    fvMsh_(fvMsh),
    solveDir_(solveDir),
    N_(fvMsh.msh().cast<rectilinearMesh>().N()),
    Nd_(Nd),
    Sd_(Sd),
    BC_(BC),
    cellSizes_(fvMsh.msh().cast<rectilinearMesh>().cellSizes()),
    D_(Nd, Zero),
    DU_(1.0 / sqr(cellSizes_[solveDir_])),
    DL_(1.0 / sqr(cellSizes_[solveDir_]))
{
    computeEigenvalues();
    computeDiagonals();
}

tridiagonalSolver::~tridiagonalSolver()
{}

void tridiagonalSolver::solve
(
    scalarBlock& f,
    scalarBlock& p,
    label start
)
{
    labelVector N(p.shape());

    // Dimension of solved direction
    label Nsolve = 0;

    // Dimensions of unsolved directions
    label N1 = 0;
    label N2 = 0;

    if (solveDir_ == 0)
    {
        Nsolve += N.x();
        N1 += N.z();
        N2 += N.y();
    }
    else if (solveDir_ == 1)
    {
        Nsolve += N.y();
        N1 += N.x();
        N2 += N.z();
    }
    else if (solveDir_ == 2)
    {
        Nsolve += N.z();
        N1 += N.x();
        N2 += N.y();
    }

    label cursor = 0;

    for (int i = 0; i < N1; i++)
    {
        for (int j = 0; j < N2; j++)
        {
            // Copy D and f so as not to overwrite them
            scalarList Dh(Nsolve, Zero);
            scalarList fh(Nsolve, Zero);

            for (int k = 0; k < Nsolve; k++)
            {
                Dh[k] = D_(cursor + k);
                fh[k] = f(cursor + k);
            }

            // Forward substitution
            for (int k = 1 + start; k < Nsolve; k++)
            {
                scalar m = DL_[k] / Dh[k-1];
                Dh[k] -= m *  DU_[k-1];
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
    scalarBlock& f,
    scalarBlock& p
)
{
    // Dimension of solved direction
    label Nsolve = 0;

    // Dimensions of unsolved directions
    label Nd1 = 0;
    label Nd2 = 0;

    if (solveDir_ == 0)
    {
        Nsolve += Nd_.x();
        Nd1 += Nd_.z();
        Nd2 += Nd_.y();
    }
    else if (solveDir_ == 1)
    {
        Nsolve += Nd_.y();
        Nd1 += Nd_.x();
        Nd2 += Nd_.z();
    }
    else if (solveDir_ == 2)
    {
        Nsolve += Nd_.z();
        Nd1 += Nd_.x();
        Nd2 += Nd_.y();
    }

    // Solve first auxiliary system
    scalarBlock u(Nd_, Zero);

    solve(f, u, 1);

    // Solve second auxiliary system
    scalarBlock v(Nd_, Zero);
    scalarBlock vf(Nd_, Zero);

    label cursor = 0;

    for (int i = 0; i < Nd1; i++)
    {
        for (int j = 0; j < Nd2; j++)
        {
            vf(cursor + 1) = -DL_[1];
            vf(cursor + Nsolve-1) = -DU_[Nsolve-1];

            cursor += Nsolve;
        }
    }

    solve(vf, v, 1);

    cursor = 0;

    // Reconstruct solution
    for (int i = 0; i < Nd1; i++)
    {
        for (int j = 0; j < Nd2; j++)
        {
            p(cursor + 0) = (f(cursor + 0) - DL_[0] * u(cursor + Nsolve-1) - DU_[0] * u(cursor + 1))
                     / (D_(cursor + 0) + DL_[0] * v(cursor + Nsolve-1) + DU_[0] * v(cursor + 1));
            for (int k = 1; k < Nsolve; k++)
            {
                p(cursor + k) = u(cursor + k) + p(cursor + 0) * v(cursor + k);
            }

            cursor += Nsolve;
        }
    }
}

void tridiagonalSolver::solve(scalarBlock& xyzPencil)
{
    scalarBlock f(xyzPencil);

    if (BC_[solveDir_] != 5)
    {
        solve(f, xyzPencil);
    }
    else
    {
        solveCyclic(f, xyzPencil);
    }
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam
