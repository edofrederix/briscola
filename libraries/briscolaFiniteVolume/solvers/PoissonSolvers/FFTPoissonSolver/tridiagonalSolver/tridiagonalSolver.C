#include "tridiagonalSolver.H"
#include "mathematicalConstants.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructors

tridiagonalSolver::tridiagonalSolver
(
    const fvMesh& fvMsh,
    labelVector Nz,
    labelVector sz,
    labelVector BC
)
:
    fvMsh_(fvMsh),
    N_(fvMsh.msh().cast<rectilinearMesh>().N()),
    Nz_(Nz),
    sz_(sz),
    BC_(BC),
    cellSizes_ (fvMsh.msh().cast<rectilinearMesh>().cellSizes())
{
    computeEigenvalues();
}

// Destructor

tridiagonalSolver::~tridiagonalSolver()
{}

// Public

void tridiagonalSolver::solve
(
    label n,
    double DU[],
    double D[],
    double DL[],
    double f[],
    double p[]
)
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

void tridiagonalSolver::solveCyclic
(
    label n,
    double DU[],
    double D[],
    double DL[],
    double f[],
    double p[]
)
{
    List<scalar> u(n, Zero);
    solve((n-1), &DU[1], &D[1], &DL[1], &f[1], &u[1]);

    List<scalar> v(n, Zero);
    List<scalar> vf(n, Zero);

    vf[1] = -DL[1];
    vf[n-1] = -DU[n-1];

    solve((n-1), &DU[1], &D[1], &DL[1], &vf[1], &v[1]);

    p[0] = (f[0] - DL[0] * u[n-1] - DU[0] * u[1])
         / (D[0] + DL[0] * v[n-1] + DU[0] * v[1]);

    for (int k = 1; k < n; k++)
    {
        p[k] = u[k] + p[0]*v[k];
    }
}

void tridiagonalSolver::solve(scalarBlock& zPencil)
{
    scalarBlock f(zPencil);

    scalar dx2 = sqr(cellSizes_[0][0]);
    scalar dy2 = sqr(cellSizes_[1][0]);
    scalarList dz2 = sqr(cellSizes_[2]);

    scalarList DU = 1.0/dz2;
    scalarList D(N_.z());
    scalarList DL = 1.0/dz2;

    switch ( BC_.z() )
    {
        case 1:
            for (int i = 0; i < Nz_.x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz_.y(); j++) // solve column by column
                {
                    D[0]        = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 3.0/dz2[0];
                    D[N_.z()-1] = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 3.0/dz2[N_.z()-1];

                    for (int k = 1; k < N_.z()-1; k++ )
                    {
                        D[k]    = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 2.0/dz2[k];
                    }
                    this->solve
                    (
                        N_.z(),
                        reinterpret_cast<double*>(DU.begin()),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(DL.begin()),
                        reinterpret_cast<double*>(&f(i,j,0)),
                        reinterpret_cast<double*>(&zPencil(i,j,0))
                    );
                }
            }
            break;

        case 2:
            for (int i = 0; i < Nz_.x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz_.y(); j++) // solve column by column
                {
                    D[0]        = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 1.0/dz2[0];
                    D[N_.z()-1] = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 1.0/dz2[N_.z()-1];
                    for (int k = 1; k < N_.z()-1; k++ )
                    {
                        D[k]    = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 2.0/dz2[k];
                    }
                    // If coefficient matrix is singular
                    if ( lambdaX_[i] == 0 && lambdaY_[j] == 0 )
                    {
                        D[0] -= 1e-10;
                    }
                    this->solve
                    (
                        N_.z(),
                        reinterpret_cast<double*>(DU.begin()),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(DL.begin()),
                        reinterpret_cast<double*>(&f(i,j,0)),
                        reinterpret_cast<double*>(&zPencil(i,j,0))
                    );
                }
            }
            break;

        case 3:
            for (int i = 0; i < Nz_.x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz_.y(); j++) // solve column by column
                {
                    D[0]        = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 3.0/dz2[0];
                    D[N_.z()-1] = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 1.0/dz2[N_.z()-1];
                    for (int k = 1; k < N_.z()-1; k++ )
                    {
                        D[k]    = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 2.0/dz2[k];
                    }
                    this->solve
                    (
                        N_.z(),
                        reinterpret_cast<double*>(DU.begin()),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(DL.begin()),
                        reinterpret_cast<double*>(&f(i,j,0)),
                        reinterpret_cast<double*>(&zPencil(i,j,0))
                    );
                }
            }
            break;

        case 4:
            for (int i = 0; i < Nz_.x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz_.y(); j++) // solve column by column
                {
                    D[0]        = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 1.0/dz2[0];
                    D[N_.z()-1] = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 3.0/dz2[N_.z()-1];
                    for (int k = 1; k < N_.z()-1; k++ )
                    {
                        D[k]    = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 2.0/dz2[k];
                    }
                    this->solve
                    (
                        N_.z(),
                        reinterpret_cast<double*>(DU.begin()),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(DL.begin()),
                        reinterpret_cast<double*>(&f(i,j,0)),
                        reinterpret_cast<double*>(&zPencil(i,j,0))
                    );
                }
            }
            break;

        case 5:
            for (int i = 0; i < Nz_.x(); i++) // loop and solve row by row
            {
                for (int j = 0; j < Nz_.y(); j++) // solve column by column
                {
                    for (int k = 0; k < N_.z(); k++ )
                    {
                        D[k] = lambdaX_[i]/dx2 + lambdaY_[j]/dy2 - 2.0/dz2[k];
                    }
                    // If coefficient matrix is singular
                    if ( lambdaX_[i] == 0 && lambdaY_[j] == 0 )
                    {
                        D[0] -= 1e-10;
                    }
                    this->solveCyclic
                    (
                        N_.z(),
                        reinterpret_cast<double*>(DU.begin()),
                        reinterpret_cast<double*>(D.begin()),
                        reinterpret_cast<double*>(DL.begin()),
                        reinterpret_cast<double*>(&f(i,j,0)),
                        reinterpret_cast<double*>(&zPencil(i,j,0))
                    );
                }
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }

}

void tridiagonalSolver::computeEigenvalues()
{
    using Foam::pow;
    using Foam::sin;

    const scalar pi(constant::mathematical::pi);

    lambdaX_ = scalarList(Nz_.x(), Zero);
    lambdaY_ = scalarList(Nz_.y(), Zero);

    switch ( BC_.x() )
    {
        case 1:
            for ( int i = 0; i < Nz_.x(); i++ )
            {
                lambdaX_[i] = - 4.0
                * pow
                (
                    sin( (sz_.x()+i+1) * pi / (2.0 * N_.x()) ),
                    2.0
                );
            }
            break;

        case 2:
            for ( int i = 0; i < Nz_.x(); i++ )
            {
                lambdaX_[i] = - 4.0
                * pow
                (
                    sin( (sz_.x()+i) * pi / (2.0 * N_.x()) ),
                    2.0
                );
            }
            break;

        case 3:
        case 4:
            for ( int i = 0; i < Nz_.x(); i++ )
            {
                lambdaX_[i] = - 4.0
                * pow
                (
                    sin( (2.0*(sz_.x()+i+1)-1.0) * pi / (4.0 * N_.x()) ),
                    2.0
                );
            }
            break;

        case 5:
            for ( int i = 0; i < Nz_.x(); i++ )
            {
                lambdaX_[i] = -4.0
                * pow
                (
                    sin( (sz_.x()+i) * pi / N_.x() ),
                    2.0
                );
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }


    switch ( BC_.y() )
    {
        case 1:
            for ( int j = 0; j < Nz_.y(); j++ )
            {
                lambdaY_[j] = - 4.0
                * pow
                (
                    sin( (sz_.y()+j+1) * pi / (2.0 * N_.y()) ),
                    2.0
                );
            }
            break;

        case 2:
            for ( int j = 0; j < Nz_.y(); j++ )
            {
                lambdaY_[j] = - 4.0
                * pow
                (
                    sin( (sz_.y()+j) * pi / (2.0 * N_.y()) ),
                    2.0
                );
            }
            break;

        case 3:
        case 4:
            for ( int j = 0; j < Nz_.y(); j++ )
            {
                lambdaY_[j] = - 4.0
                * pow
                (
                    sin( (2.0*(sz_.y()+j+1)-1.0) * pi / (4.0 * N_.y()) ),
                    2.0
                );
            }
            break;

        case 5:
            for ( int j = 0; j < Nz_.y(); j++ )
            {
                lambdaY_[j] = -4.0
                * pow
                (
                    sin( (sz_.y()+j) * pi / N_.y() ),
                    2.0
                );
            }
            break;

        default:
            FatalError << " Unknown boundary condition "
            << endl << abort(FatalError);
    }
}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam
