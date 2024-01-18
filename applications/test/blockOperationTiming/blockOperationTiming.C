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

    fvMesh fvMsh(meshDict, runTime);

    const labelVector N = fvMsh.msh().bricks()[0].N();

    // Internal block creation

    vectorBlock b1(N);
    vectorBlock b2(N);

    forAllBlock(b1, i, j, k)
    {
        b1(i,j,k) = vector(i,j,k);
        b2(i,j,k) = vector(i*2,j*3,-k*4);
    }

    tensorBlock b3(b2.shape());

    // Internal block indexed loop

    auto t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllBlock(b3, i, j, k)
            b3(i,j,k) = b1(i,j,k) * b2(i,j,k);

    auto t2 = high_resolution_clock::now();

    const label tInternalBlockIndexed =
        duration_cast<milliseconds>(t2 - t1).count();

    Info<< "Internal block indexed loop =                       "
        << tInternalBlockIndexed
        << " ms" << endl;

    // Internal block linear loop

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllBlockLinear(b3, i)
            b3(i) = b1(i) * b2(i);

    t2 = high_resolution_clock::now();

    const label tInternalBlockLinear =
        duration_cast<milliseconds>(t2 - t1).count();

    Info<< "Internal block linear loop =                        "
        << tInternalBlockLinear
        << " ms" << endl;

    const scalar costOfIndexing =
        scalar(tInternalBlockIndexed)/scalar(tInternalBlockLinear);

    // Padded block creation

    vectorBlock p1(N+P*2*unitXYZ);
    vectorBlock p2(N+P*2*unitXYZ);

    for (label i = P; i < p1.l()-P; i++)
    for (label j = P; j < p1.m()-P; j++)
    for (label k = P; k < p1.n()-P; k++)
    {
        p1(i,j,k) = vector(i-P,j-P,k-P);
        p2(i,j,k) = vector((i-P)*2,(j-P)*3,-(k-P)*4);
    }

    tensorBlock p3(p2.shape());

    // Padded block indexed loop over internal cells

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        for (label i = P; i < p1.l()-P; i++)
            for (label j = P; j < p1.m()-P; j++)
                for (label k = P; k < p1.n()-P; k++)
                    p3(i,j,k) = p1(i,j,k) * p2(i,j,k);

    t2 = high_resolution_clock::now();

    const label tPaddedBlockIndexed =
        duration_cast<milliseconds>(t2 - t1).count();

    Info<< "Padded block indexed loop over internal cells =     "
        << tPaddedBlockIndexed
        << " ms" << endl;

    const scalar costOfPadding =
        scalar(tPaddedBlockIndexed)/scalar(tInternalBlockIndexed);

    // Full padded block linear loop

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllBlockLinear(p3, i)
            p3(i) = p1(i) * p2(i);

    t2 = high_resolution_clock::now();

    {
        const scalar t = duration_cast<milliseconds>(t2 - t1).count();
        const label expectation =
            scalar(p1.size())/scalar(b1.size())*tInternalBlockLinear;

        Info<< "Full padded block linear loop =                     "
            << t << " ms\t"
            << "expectation = " << expectation << " ms\t"
            << "relative = " << round(t/expectation*100) << " %" << endl;
    }

    // Mesh field creation

    colocatedVectorField f1("f1", fvMsh);
    colocatedVectorField f2("f2", fvMsh);

    forAllCells(f1, i, j, k)
    {
        f1(i,j,k) = vector(i,j,k);
        f2(i,j,k) = vector(i*2,j*3,-k*4);
    }

    colocatedTensorField f3("f3", fvMsh);

    // Internal cell indexed loop referencing fields

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllCells(f3, i, j, k)
            f3(i,j,k) = f1(i,j,k) * f2(i,j,k);

    t2 = high_resolution_clock::now();

    {
        const scalar t = duration_cast<milliseconds>(t2 - t1).count();
        const label expectation = tInternalBlockIndexed*costOfPadding;

        Info<< "Internal cell indexed loop referencing fields =     "
            << t << " ms\t"
            << "expectation = " << expectation << " ms\t"
            << "relative = " << round(t/expectation*100) << " %" << endl;
    }

    // Internal cell indexed loop referencing levels

    const colocatedVectorLevel& lf1 = f1[0];
    const colocatedVectorLevel& lf2 = f2[0];
    colocatedTensorLevel& lf3 = f3[0];

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllCells(lf3, i, j, k)
            lf3(i,j,k) = lf1(i,j,k) * lf2(i,j,k);

    t2 = high_resolution_clock::now();

    {
        const scalar t = duration_cast<milliseconds>(t2 - t1).count();
        const label expectation = tInternalBlockIndexed*costOfPadding;

        Info<< "Internal cell indexed loop referencing levels =     "
            << t << " ms\t"
            << "expectation = " << expectation << " ms\t"
            << "relative = " << round(t/expectation*100) << " %" << endl;
    }

    // Internal cell indexed loop referencing directions

    const colocatedVectorDirection& df1 = f1.direction();
    const colocatedVectorDirection& df2 = f2.direction();
    colocatedTensorDirection& df3 = f3.direction();

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllCells(df3, i, j, k)
            df3(i,j,k) = df1(i,j,k) * df2(i,j,k);

    t2 = high_resolution_clock::now();

    {
        const scalar t = duration_cast<milliseconds>(t2 - t1).count();
        const label expectation = tInternalBlockIndexed*costOfPadding;

        Info<< "Internal cell indexed loop referencing directions = "
            << t << " ms\t"
            << "expectation = " << expectation << " ms\t"
            << "relative = " << round(t/expectation*100) << " %" << endl;
    }

    // Full block indexed loop

    const vectorBlock& bf1 = f1.B();
    const vectorBlock& bf2 = f2.B();
    tensorBlock& bf3 = f3.B();

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllBlock(bf1, i, j, k)
            bf3(i,j,k) = bf1(i,j,k) * bf2(i,j,k);

    t2 = high_resolution_clock::now();

    {
        const scalar t = duration_cast<milliseconds>(t2 - t1).count();
        const label expectation =
            scalar(p1.size())/scalar(b1.size())
          * tInternalBlockLinear*costOfIndexing;

        Info<< "Full block indexed loop =                           "
            << t << " ms\t"
            << "expectation = " << expectation << " ms\t"
            << "relative = " << round(t/expectation*100) << " %" << endl;
    }

    // Full block linear loop

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        forAllBlockLinear(bf3, i)
            bf3(i) = bf1(i) * bf2(i);

    t2 = high_resolution_clock::now();

    {
        const scalar t = duration_cast<milliseconds>(t2 - t1).count();
        const label expectation =
            scalar(p1.size())/scalar(b1.size())*tInternalBlockLinear;

        Info<< "Full block linear loop =                            "
            << t << " ms\t"
            << "expectation = " << expectation << " ms\t"
            << "relative = " << round(t/expectation*100) << " %" << endl;
    }

    // Full block operation with tmp

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        f3.B() = f1.B() * f2.B();

    t2 = high_resolution_clock::now();

    {
        const scalar t = duration_cast<milliseconds>(t2 - t1).count();
        const label expectation =
            scalar(p1.size())/scalar(b1.size())*tInternalBlockLinear;

        Info<< "Full block operation with tmp (slow) =              "
            << t << " ms\t"
            << "expectation = " << expectation << " ms\t"
            << "relative = " << round(t/expectation*100) << " %" << endl;
    }

    // Full block operation without tmp

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        outer(f3.B(), f1.B(), f2.B());

    t2 = high_resolution_clock::now();

    {
        const scalar t = duration_cast<milliseconds>(t2 - t1).count();
        const label expectation =
            scalar(p1.size())/scalar(b1.size())*tInternalBlockLinear;

        Info<< "Full block operation without tmp =                  "
            << t << " ms\t"
            << "expectation = " << expectation << " ms\t"
            << "relative = " << round(t/expectation*100) << " %" << endl;
    }

    // Full direction operation with tmp

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        f3.direction() = f1.direction() * f2.direction();

    t2 = high_resolution_clock::now();

    {
        const scalar t = duration_cast<milliseconds>(t2 - t1).count();
        const label expectation =
            scalar(p1.size())/scalar(b1.size())*tInternalBlockLinear;

        Info<< "Full direction operation with tmp (slow) =          "
            << t << " ms\t"
            << "expectation = " << expectation << " ms\t"
            << "relative = " << round(t/expectation*100) << " %" << endl;
    }

    // Full direction operation without tmp

    t1 = high_resolution_clock::now();

    for (label iter = 0; iter < Niter; iter++)
        outer(f3.direction(), f1.direction(), f2.direction());

    t2 = high_resolution_clock::now();

    {
        const scalar t = duration_cast<milliseconds>(t2 - t1).count();
        const label expectation =
            scalar(p1.size())/scalar(b1.size())*tInternalBlockLinear;

        Info<< "Full direction operation without tmp =              "
            << t << " ms\t"
            << "expectation = " << expectation << " ms\t"
            << "relative = " << round(t/expectation*100) << " %" << endl;
    }
}
