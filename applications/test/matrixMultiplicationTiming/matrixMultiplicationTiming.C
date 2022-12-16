#include "arguments.H"
#include "Time.H"

#include "fileOperation.H"
#include "OSspecific.H"

#include "IFstream.H"
#include "OFstream.H"

#include "colocatedFields.H"

#include "fvMesh.H"

#include <chrono>

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
    using std::chrono::high_resolution_clock;
    using std::chrono::duration_cast;
    using std::chrono::duration;
    using std::chrono::milliseconds;

    #include "createParallelBriscolaCase.H"
    #include "createBriscolaTime.H"

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

    const label Niter = 200;

    fvMesh fvMsh(meshDict, runTime);

    // Fields, levels and directions

    {
        colocatedVectorField  F("F", fvMsh);
        colocatedVectorField  G("G", fvMsh);
        colocatedStencilField A("A", fvMsh);

        F = Zero;

        forAll(A, l)
        forAll(A[l], d)
        forAllBlock(A[l][d], i, j, k)
        {
            F[l][d](i,j,k) = vector(i,j,k);
            A[l][d](i,j,k) = stencil(-6, 1, 1, 1, 1, 1, 1);
        }

        // Matrix multiplication for a whole field. Will be slower because it
        // multiplies all levels and directions

        auto t1 = high_resolution_clock::now();

        for (label iter = 0; iter < Niter; iter++)
            Amul(G, A, F);

        auto t2 = high_resolution_clock::now();

        Info<< "Amul field time = "
            << duration_cast<milliseconds>(t2 - t1).count()
            << " ms" << endl;

        // Matrix multiplication for a level. Should take about the same time as for
        // a direction, because the level only has one direction (colocated)

        t1 = high_resolution_clock::now();

        for (label iter = 0; iter < Niter; iter++)
            Amul(G[0], A[0], F[0]);

        t2 = high_resolution_clock::now();

        Info<< "Amul level time = "
            << duration_cast<milliseconds>(t2 - t1).count()
            << " ms" << endl;

        // Matrix multiplcation for a direction

        t1 = high_resolution_clock::now();

        for (label iter = 0; iter < Niter; iter++)
            Amul(G[0][0], A[0][0], F[0][0]);

        t2 = high_resolution_clock::now();

        Info<< "Amul direction time = "
            << duration_cast<milliseconds>(t2 - t1).count()
            << " ms" << endl;
    }

    // Blocks

    {
        vectorBlock bF(fvMsh.N()+unitXYZ*2);
        vectorBlock bG(fvMsh.N()+unitXYZ*2);
        stencilBlock bA(fvMsh.N()+unitXYZ*2);

        forAllBlock(bF, i, j, k)
        {
            bF(i,j,k) = vector(i,j,k);
            bA(i,j,k) = stencil(-6, 1, 1, 1, 1, 1, 1);
        }

        // Custom matrix multiplcation for a block with the same size as a
        // direction. Should be slightly faster than matrix multiplcation for a
        // direction, because the ghost cells offset is not computed for every
        // call to operator()

        auto t1 = high_resolution_clock::now();

        for (label iter = 0; iter < Niter; iter++)
        for (label i = 1; i < bF.l()-1; i++)
        for (label j = 1; j < bF.m()-1; j++)
        for (label k = 1; k < bF.n()-1; k++)
        {
            bG(i,j,k) =
                bA(i,j,k).center() * bF(i,  j,  k  )
              + bA(i,j,k).left()   * bF(i-1,j,  k  )
              + bA(i,j,k).right()  * bF(i+1,j,  k  )
              + bA(i,j,k).bottom() * bF(i,  j-1,k  )
              + bA(i,j,k).top()    * bF(i,  j+1,k  )
              + bA(i,j,k).aft()    * bF(i,  j,  k-1)
              + bA(i,j,k).fore()   * bF(i,  j,  k+1);
        }

        auto t2 = high_resolution_clock::now();

        Info<< "Block multiplication = "
            << duration_cast<milliseconds>(t2 - t1).count()
            << " ms" << endl;

        scalarBlock bc(bF.shape(), -6);
        scalarBlock bl(bF.shape(), 1);
        scalarBlock br(bF.shape(), 1);
        scalarBlock bb(bF.shape(), 1);
        scalarBlock bt(bF.shape(), 1);
        scalarBlock ba(bF.shape(), 1);
        scalarBlock bf(bF.shape(), 1);

        // Matrix multiplication by split matrix coefficient fields

        t1 = high_resolution_clock::now();

        for (label iter = 0; iter < Niter; iter++)
        for (label i = 1; i < bF.l()-1; i++)
        for (label j = 1; j < bF.m()-1; j++)
        for (label k = 1; k < bF.n()-1; k++)
        {
            bG(i,j,k) =
                bc(i,j,k) * bF(i,  j,  k  )
              + bl(i,j,k) * bF(i-1,j,  k  )
              + br(i,j,k) * bF(i+1,j,  k  )
              + bb(i,j,k) * bF(i,  j-1,k  )
              + bt(i,j,k) * bF(i,  j+1,k  )
              + ba(i,j,k) * bF(i,  j,  k-1)
              + bf(i,j,k) * bF(i,  j,  k+1);
        }

        t2 = high_resolution_clock::now();

        Info<< "Block multiplication with split matrix = "
            << duration_cast<milliseconds>(t2 - t1).count()
            << " ms" << endl;
    }


    // Low level c++ implementation

    {
        const int Nx = fvMsh.N().x() + 2;
        const int Ny = fvMsh.N().y() + 2;
        const int Nz = fvMsh.N().z() + 2;
        const int Nt = Nx*Ny*Nz;

        auto F = new double[Nt][3];
        auto G = new double[Nt][3];
        auto A = new double[Nt][6];

        auto c = new double[Nt];
        auto l = new double[Nt];
        auto r = new double[Nt];
        auto b = new double[Nt];
        auto t = new double[Nt];
        auto a = new double[Nt];
        auto f = new double[Nt];

        for (int i = 1; i < Nx-1; i++)
        for (int j = 1; j < Ny-1; j++)
        for (int k = 1; k < Nz-1; k++)
        {
            int ijk = i*Ny*Nz + j*Nz + k;

            F[ijk][0] = i;
            F[ijk][1] = j;
            F[ijk][2] = k;

            c[ijk] = -6;
            l[ijk] = 1;
            r[ijk] = 1;
            b[ijk] = 1;
            t[ijk] = 1;
            a[ijk] = 1;
            f[ijk] = 1;

            A[ijk][0] = c[ijk];
            A[ijk][0] = l[ijk];
            A[ijk][0] = r[ijk];
            A[ijk][0] = b[ijk];
            A[ijk][0] = t[ijk];
            A[ijk][0] = a[ijk];
            A[ijk][0] = f[ijk];
        }

        auto t1 = high_resolution_clock::now();

        for (int iter = 0; iter < Niter; iter++)
        for (int i = 1; i < Nx-1; i++)
        for (int j = 1; j < Ny-1; j++)
        for (int k = 1; k < Nz-1; k++)
        {
            int ijk = i*Ny*Nz + j*Nz + k;

            for (int d = 0; d < 3; d++)
                G[ijk][d] =
                    A[ijk][0] * F[(i  )*Ny*Nz + (j  )*Nz + (k  )][d]
                  + A[ijk][1] * F[(i-1)*Ny*Nz + (j  )*Nz + (k  )][d]
                  + A[ijk][2] * F[(i+1)*Ny*Nz + (j  )*Nz + (k  )][d]
                  + A[ijk][3] * F[(i  )*Ny*Nz + (j-1)*Nz + (k  )][d]
                  + A[ijk][4] * F[(i  )*Ny*Nz + (j+1)*Nz + (k  )][d]
                  + A[ijk][5] * F[(i  )*Ny*Nz + (j  )*Nz + (k-1)][d]
                  + A[ijk][6] * F[(i  )*Ny*Nz + (j  )*Nz + (k+1)][d];
        }

        auto t2 = high_resolution_clock::now();

        Info<< "C++ array multiplication = "
            << duration_cast<milliseconds>(t2 - t1).count()
            << " ms" << endl;

        t1 = high_resolution_clock::now();

        for (int iter = 0; iter < Niter; iter++)
        for (int i = 1; i < Nx-1; i++)
        for (int j = 1; j < Ny-1; j++)
        for (int k = 1; k < Nz-1; k++)
        {
            int ijk = i*Ny*Nz + j*Nz + k;

            for (int d = 0; d < 3; d++)
                G[ijk][d] =
                    c[ijk] * F[(i  )*Ny*Nz + (j  )*Nz + (k  )][d]
                  + l[ijk] * F[(i-1)*Ny*Nz + (j  )*Nz + (k  )][d]
                  + r[ijk] * F[(i+1)*Ny*Nz + (j  )*Nz + (k  )][d]
                  + b[ijk] * F[(i  )*Ny*Nz + (j-1)*Nz + (k  )][d]
                  + t[ijk] * F[(i  )*Ny*Nz + (j+1)*Nz + (k  )][d]
                  + a[ijk] * F[(i  )*Ny*Nz + (j  )*Nz + (k-1)][d]
                  + f[ijk] * F[(i  )*Ny*Nz + (j  )*Nz + (k+1)][d];
        }

        t2 = high_resolution_clock::now();

        Info<< "C++ array multiplication with split matrix = "
            << duration_cast<milliseconds>(t2 - t1).count()
            << " ms" << endl;
    }
}
