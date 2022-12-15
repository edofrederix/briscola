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

    const label Niter = 1000;

    fvMesh fvMsh(meshDict, runTime);

    colocatedVectorField F("F", fvMsh.metrics<colocated>().cellCenters());
    colocatedVectorField G("G", fvMsh);
    colocatedStencilField A("A", fvMsh);

    forAll(A, l)
        forAll(A[l], d)
            forAllBlock(A[l][d], i, j, k)
                A[l][d](i,j,k) = stencil(-6, 1, 1, 1, 1, 1, 1);

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


    vectorBlock bF(A[0][0].B().shape());
    vectorBlock bG(A[0][0].B().shape());
    stencilBlock bA(A[0][0].B().shape());

    forAllBlock(bF, i, j, k)
    {
        bF(i,j,k) = F[0][0](i,j,k);
        bA(i,j,k) = stencil(-6, 1, 1, 1, 1, 1, 1);
    }

    // Custom matrix multiplcation for a block with the same size as a
    // direction. Should be slightly faster than matrix multiplcation for a
    // direction, because the ghost cells offset is not computed for every call
    // to operator()

    t1 = high_resolution_clock::now();

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

    t2 = high_resolution_clock::now();

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
