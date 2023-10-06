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

    const label Niter = 100;
    const label P = 1;
    const label N = 64;

    fvMesh fvMsh(meshDict, runTime);


    // Normal block creation

    vectorBlock b1(N, N*2, N*3);
    vectorBlock b2(N, N*2, N*3);

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = vector(i,j,k);
        b2(i,j,k) = vector(i*2,j*3,-k*4);
    }

    tensorBlock b3(b2.shape());

    // Normal block calculation: loop over elements in block in a structured
    // manner

    auto t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllBlock(b3, i, j, k)
            b3(i,j,k) = b1(i,j,k) * b2(i,j,k);

    auto t2 = high_resolution_clock::now();

    Info<< "Normal block calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Normal block linear calculation: loop over elements in a block in a
    // linear manner. This should be slightly cheaper than the structured
    // method.

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllBlockLinear(b3, i)
            b3(i) = b1(i) * b2(i);

    t2 = high_resolution_clock::now();

    Info<< "Normal block linear calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Padded block creation

    vectorBlock p1(N+P*2, N*2+P*2, N*3+P*2);
    vectorBlock p2(N+P*2, N*2+P*2, N*3+P*2);

    for (label i = P; i < p1.l()-P; i++)
    for (label j = P; j < p1.m()-P; j++)
    for (label k = P; k < p1.n()-P; k++)
    {
        p1(i,j,k) = vector(i-P,j-P,k-P);
        p2(i,j,k) = vector((i-P)*2,(j-P)*3,-(k-P)*4);
    }

    tensorBlock p3(p2.shape());

    // Padded block calculation: account for padding in block

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        for (label i = P; i < p1.l()-P; i++)
            for (label j = P; j < p1.m()-P; j++)
                for (label k = P; k < p1.n()-P; k++)
                    p3(i,j,k) = p1(i,j,k) * p2(i,j,k);

    t2 = high_resolution_clock::now();

    Info<< "Padded block calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Padded block linear calculation: don't account for padding but perform
    // linear operation on the full block instead.

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllBlockLinear(p3, i)
            p3(i) = p1(i) * p2(i);

    t2 = high_resolution_clock::now();

    Info<< "Padded block linear calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Mesh field creation

    colocatedVectorField f1("f1", fvMsh);
    colocatedVectorField f2("f2", fvMsh);

    forAllCells(f1, i, j, k)
    {
        f1(i,j,k) = vector(i,j,k);
        f2(i,j,k) = vector(i*2,j*3,-k*4);
    }

    colocatedTensorField f3("f3", fvMsh);

    // Field reference looped: use the level and direction access operators,
    // which gives a slight overhead.

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllCells(f3, i, j, k)
            f3(i,j,k) = f1(i,j,k) * f2(i,j,k);

    t2 = high_resolution_clock::now();

    Info<< "Field reference looped calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Level reference looped: use the direction access operator, which gives a
    // slight overhead.

    const colocatedVectorLevel& lf1 = f1[0];
    const colocatedVectorLevel& lf2 = f2[0];
    colocatedTensorLevel& lf3 = f3[0];

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllCells(lf3, i, j, k)
            lf3(i,j,k) = lf1(i,j,k) * lf2(i,j,k);

    t2 = high_resolution_clock::now();

    Info<< "Field level reference looped calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Direction reference looped: avoid level/direction access operators and
    // access elements in the direction directly.

    const colocatedVectorDirection& df1 = f1.direction();
    const colocatedVectorDirection& df2 = f2.direction();
    colocatedTensorDirection& df3 = f3.direction();

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllCells(df3, i, j, k)
            df3(i,j,k) = df1(i,j,k) * df2(i,j,k);

    t2 = high_resolution_clock::now();

    Info<< "Field direction reference looped calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Block reference looped: perform structured loop over block data directly.
    // This gives a bit of overhead because the ghost cells are also taken into
    // account.

    const vectorBlock& bf1 = f1.B();
    const vectorBlock& bf2 = f2.B();
    tensorBlock& bf3 = f3.B();

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllBlock(bf3, i, j, k)
            bf3(i,j,k) = bf1(i,j,k) * bf2(i,j,k);

    t2 = high_resolution_clock::now();

    Info<< "Field block reference looped calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Block reference linear looped: linear loop over block data directly.
    // Gives a bit of overhead because ghost cells are also taken into account,
    // but should also be faster because of the linear traversal.

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllBlockLinear(bf3, i)
            bf3(i) = bf1(i) * bf2(i);

    t2 = high_resolution_clock::now();

    Info<< "Field block reference linear looped calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Block: this is quite expensive because at each operation a new temporary
    // block is created. Avoid if possible.

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        f3.B() = f1.B() * f2.B();

    t2 = high_resolution_clock::now();

    Info<< "Field block operator calculation time (expensive) = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Block function: use the equivalent of the operator, which avoids creation
    // of a temporary block and is thus faster.

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        outer(f3.B(), f1.B(), f2.B());

    t2 = high_resolution_clock::now();

    Info<< "Field block function calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Direct operator: perform field operator, which is expensive because
    // temporary memory is created.

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        f3.direction() = f1.direction() * f2.direction();

    t2 = high_resolution_clock::now();

    Info<< "Field direct operator calculation time (expensive)= "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;

    // Direct function: use the equivalent of the operator, which avoids
    // creation of a temporary direction and is thus faster.

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        outer(f3.direction(), f1.direction(), f2.direction());

    t2 = high_resolution_clock::now();

    Info<< "Field direct function calculation time = "
        << duration_cast<milliseconds>(t2 - t1).count()
        << " ms" << endl;
}
