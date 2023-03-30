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

    lambdaX_ = scalarList(Nz_.x(), Zero);
    lambdaY_ = scalarList(Nz_.y(), Zero);

    switch ( BC_.x() )
    {
        case 1:
            for (int i = 0; i < Nz_.x(); i++)
            {
                lambdaX_[i] = - 4.0
                * sqr(sin((Sz_.x()+i+1) * pi / (2.0 * N_.x())));
            }
            break;

        case 2:
            for (int i = 0; i < Nz_.x(); i++)
            {
                lambdaX_[i] = - 4.0
                * sqr(sin((Sz_.x()+i) * pi / (2.0 * N_.x())));
            }
            break;

        case 3:
        case 4:
            for (int i = 0; i < Nz_.x(); i++)
            {
                lambdaX_[i] = - 4.0
                * sqr(sin((2.0*(Sz_.x()+i+1)-1.0) * pi / (4.0 * N_.x())));
            }
            break;

        case 5:
            for (int i = 0; i < Nz_.x(); i++)
            {
                lambdaX_[i] = -4.0
                * sqr(sin((Sz_.x()+i) * pi / N_.x()));
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }

    switch ( BC_.y() )
    {
        case 1:
            for (int j = 0; j < Nz_.y(); j++)
            {
                lambdaY_[j] = - 4.0
                * sqr(sin((Sz_.y()+j+1) * pi / (2.0 * N_.y())));
            }
            break;

        case 2:
            for (int j = 0; j < Nz_.y(); j++)
            {
                lambdaY_[j] = - 4.0
                * sqr(sin((Sz_.y()+j) * pi / (2.0 * N_.y())));
            }
            break;

        case 3:
        case 4:
            for (int j = 0; j < Nz_.y(); j++)
            {
                lambdaY_[j] = - 4.0
                * sqr(sin((2.0*(Sz_.y()+j+1)-1.0) * pi / (4.0 * N_.y())));
            }
            break;

        case 5:
            for (int j = 0; j < Nz_.y(); j++)
            {
                lambdaY_[j] = -4.0
                * sqr(sin((Sz_.y()+j) * pi / N_.y()));
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }
}

void tridiagonalSolver::computeDiagonals()
{
    scalar dx2 = sqr(cellSizes_[0][0]);
    scalar dy2 = sqr(cellSizes_[1][0]);
    scalarList dz2 = sqr(cellSizes_[2]);

    for (int i = 0; i < Nz_.x(); i++)
    {
        for (int j = 0; j < Nz_.y(); j++)
        {
            switch ( BC_.z() )
            {
                case 1:
                    D_(i,j,0) = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 3.0/dz2[0];
                    D_(i,j,N_.z()-1) = lambdaX_[i]/dx2
                                     + lambdaY_[j]/dy2 - 3.0/dz2[N_.z()-1];
                    break;

                case 2:
                    D_(i,j,0) = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 1.0/dz2[0];
                    D_(i,j,N_.z()-1) = lambdaX_[i]/dx2
                                     + lambdaY_[j]/dy2 - 1.0/dz2[N_.z()-1];
                    // If coefficient matrix is singular
                    if ( lambdaX_[i] == 0 && lambdaY_[j] == 0 )
                    {
                        D_(i,j,0) -= 1e-10;
                    }
                    break;

                case 3:
                    D_(i,j,0) = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 3.0/dz2[0];
                    D_(i,j,N_.z()-1) = lambdaX_[i]/dx2
                                     + lambdaY_[j]/dy2 - 1.0/dz2[N_.z()-1];
                    break;

                case 4:
                    D_(i,j,0) = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 1.0/dz2[0];
                    D_(i,j,N_.z()-1) = lambdaX_[i]/dx2
                                     + lambdaY_[j]/dy2 - 3.0/dz2[N_.z()-1];
                    break;

                case 5:
                    D_(i,j,0) = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 2.0/dz2[0];
                    D_(i,j,N_.z()-1) = lambdaX_[i]/dx2
                                     + lambdaY_[j]/dy2 - 2.0/dz2[N_.z()-1];
                    // If coefficient matrix is singular
                    if ( lambdaX_[i] == 0 && lambdaY_[j] == 0 )
                    {
                        D_(i,j,0) -= 1e-10;
                    }
                    break;
            }

            for (int k = 1; k < Nz_.z() - 1; k++)
            {
                D_(i,j,k) = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 2.0/dz2[k];
            }
        }
    }
}

tridiagonalSolver::tridiagonalSolver
(
    const fvMesh& fvMsh,
    labelVector Nz,
    labelVector Sz,
    labelVector BC
)
:
    fvMsh_(fvMsh),
    N_(fvMsh.msh().cast<rectilinearMesh>().N()),
    Nz_(Nz),
    Sz_(Sz),
    BC_(BC),
    cellSizes_(fvMsh.msh().cast<rectilinearMesh>().cellSizes()),
    D_(Nz, Zero),
    DU_(1.0 / sqr(cellSizes_[2])),
    DL_(1.0 / sqr(cellSizes_[2]))
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

    for (int i = 0; i < N.x(); i++)
    {
        for (int j = 0; j < N.y(); j++)
        {
            // Copy D and f so as not to overwrite them
            scalarList Dh(N.z(), Zero);
            scalarList fh(N.z(), Zero);

            for (int k = 0; k < N.z(); k++)
            {
                Dh[k] = D_(i,j,k);
                fh[k] = f(i,j,k);
            }

            // Forward substitution
            for (int k = 1 + start; k < N.z(); k++)
            {
                scalar m = DL_[k] / Dh[k-1];
                Dh[k] -= m *  DU_[k-1];
                fh[k] -= m * fh[k-1];
            }

            // Backward substitution
            p(i,j,N.z()-1) = fh[N.z()-1] / Dh[N.z()-1];
            for (int k = N.z()-2; k >= start; k--)
            {
                p(i,j,k) = (fh[k] - DU_[k] * p(i,j,k+1)) / Dh[k];
            }
        }
    }
}

void tridiagonalSolver::solveCyclic
(
    scalarBlock& f,
    scalarBlock& p
)
{
    // Solve first auxiliary system
    scalarBlock u(Nz_, Zero);

    solve(f, u, 1);

    // Solve second auxiliary system
    scalarBlock v(Nz_, Zero);
    scalarBlock vf(Nz_, Zero);

    for (int i = 0; i < Nz_.x(); i++)
    {
        for (int j = 0; j < Nz_.y(); j++)
        {
            vf(i,j,1) = -DL_[1];
            vf(i,j,Nz_.z()-1) = -DU_[Nz_.z()-1];
        }
    }

    solve(vf, v, 1);

    // Reconstruct solution
    for (int i = 0; i < Nz_.x(); i++)
    {
        for (int j = 0; j < Nz_.y(); j++)
        {
            p(i,j,0) = (f(i,j,0) - DL_[0] * u(i,j,N_.z()-1) - DU_[0] * u(i,j,1))
                     / (D_(i,j,0) + DL_[0] * v(i,j,N_.z()-1) + DU_[0] * v(i,j,1));
            for (int k = 1; k < Nz_.z(); k++)
            {
                p(i,j,k) = u(i,j,k) + p(i,j,0) * v(i,j,k);
            }
        }
    }
}

void tridiagonalSolver::solve(scalarBlock& zPencil)
{
    scalarBlock f(zPencil);

    if (BC_.z() != 5)
    {
        solve(f, zPencil);
    }
    else
    {
        solveCyclic(f, zPencil);
    }
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam
