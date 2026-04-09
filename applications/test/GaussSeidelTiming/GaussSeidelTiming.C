#include "arguments.H"
#include "Time.H"

#include "colocatedFields.H"
#include "staggeredFields.H"

#include "fvMesh.H"

#include "faceFields.H"
#include "stencilFields.H"
#include "diagStencilFields.H"

#include "ticToc.H"

using namespace Foam;
using namespace briscola;
using namespace fv;

int main(int argc, char *argv[])
{
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

    fvMesh fvMsh(meshDict, runTime);

    initTicToc(32)

    colocatedScalarStencilField A("A", fvMsh);

    for (int c = 1; c < A.size(); c++)
        forAllCells(A[c], i, j, k)
            A[c](i,j,k) = c-i+j+k;

    colocatedStencilField S("S", fvMsh);

    forAllCells(S, i, j, k)
        for(int c = 0; c < 7; c++)
            S(i,j,k)[c] = c-i+j+k;

    colocatedScalarField xs("xs", fvMsh);
    colocatedScalarField bs("bs", fvMsh);

    forAllBlock(xs.B(), i, j, k)
    {
        xs.B()(i,j,k) = i-j+k;
        bs.B()(i,j,k) = -2*i-j+k;
    }

    colocatedVectorField xv("xv", fvMsh);
    colocatedVectorField bv("bv", fvMsh);

    forAllBlock(xv.B(), i, j, k)
    {
        xv.B()(i,j,k) = vector::one*(i+j-k);
        bv.B()(i,j,k) = vector::one*(-2*i+j-k);
    }

    const labelVector* offsets = stencil::offsets;

    for (int iter = 0; iter < 200; iter++)
    {
        xs = Zero;
        xv = Zero;

        // RBGS

        // Scalar

            // SoA storage

                // Outer loop

                Tic(0)

                for (int color = 0; color < 2; color++)
                {
                    forAllCells(xs, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xs(i,j,k) = bs(i,j,k);

                    for (int c = 1; c < A.size(); c++)
                        forAllCells(A[c], i, j, k)
                            if ((i + j + k) % 2 == color)
                                xs(i,j,k) -=
                                    A[c](i,j,k)
                                  * xs(labelVector(i,j,k) + offsets[c]);

                    forAllCells(xs, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xs(i,j,k) /= A[0](i,j,k);
                }

                Toc(0)

                // Inner loop

                Tic(1)

                for (int color = 0; color < 2; color++)
                {
                    forAllCells(xs, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xs(i,j,k) = bs(i,j,k);

                    forAllCells(xs, i, j, k)
                        if ((i + j + k) % 2 == color)
                            for (int c = 1; c < A.size(); c++)
                                xs(i,j,k) -=
                                    A[c](i,j,k)
                                  * xs(labelVector(i,j,k) + offsets[c]);

                    forAllCells(xs, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xs(i,j,k) /= A[0](i,j,k);
                }

                Toc(1)

            // AoS storage

                // Outer loop

                Tic(2)

                for (int color = 0; color < 2; color++)
                {

                    forAllCells(xs, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xs(i,j,k) = bs(i,j,k);

                    for (int c = 1; c < A.size(); c++)
                        forAllCells(S, i, j, k)
                            if ((i + j + k) % 2 == color)
                                xs(i,j,k) -=
                                    S(i,j,k)[c]
                                  * xs(labelVector(i,j,k) + offsets[c]);

                    forAllCells(xs, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xs(i,j,k) /= S(i,j,k)[0];
                }

                Toc(2)

                // Inner loop

                Tic(3)

                for (int color = 0; color < 2; color++)
                {
                    forAllCells(xs, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xs(i,j,k) = bs(i,j,k);

                    forAllCells(xs, i, j, k)
                        if ((i + j + k) % 2 == color)
                            for (int c = 1; c < A.size(); c++)
                                xs(i,j,k) -=
                                    S(i,j,k)[c]
                                  * xs(labelVector(i,j,k) + offsets[c]);

                    forAllCells(xs, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xs(i,j,k) /= S(i,j,k)[0];
                }

                Toc(3)



        // Vector

            // SoA storage

                // Outer loop

                Tic(4)

                for (int color = 0; color < 2; color++)
                {
                    forAllCells(xv, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xv(i,j,k) = bv(i,j,k);

                    for (int c = 1; c < A.size(); c++)
                        forAllCells(A[c], i, j, k)
                            if ((i + j + k) % 2 == color)
                                xv(i,j,k) -=
                                    A[c](i,j,k)
                                  * xv(labelVector(i,j,k) + offsets[c]);

                    forAllCells(xv, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xv(i,j,k) /= A[0](i,j,k);
                }

                Toc(4)

                // Inner loop

                Tic(5)

                for (int color = 0; color < 2; color++)
                {
                    forAllCells(xv, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xv(i,j,k) = bv(i,j,k);

                    forAllCells(xv, i, j, k)
                        if ((i + j + k) % 2 == color)
                            for (int c = 1; c < A.size(); c++)
                                xv(i,j,k) -=
                                    A[c](i,j,k)
                                  * xv(labelVector(i,j,k) + offsets[c]);

                    forAllCells(xv, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xv(i,j,k) /= A[0](i,j,k);
                }

                Toc(5)

            // AoS storage

                // Outer loop

                Tic(6)

                for (int color = 0; color < 2; color++)
                {
                    forAllCells(xv, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xv(i,j,k) = bv(i,j,k);

                    for (int c = 1; c < A.size(); c++)
                        forAllCells(S, i, j, k)
                            if ((i + j + k) % 2 == color)
                                xv(i,j,k) -=
                                    S(i,j,k)[c]
                                  * xv(labelVector(i,j,k) + offsets[c]);

                    forAllCells(xv, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xv(i,j,k) /= S(i,j,k)[0];
                }

                Toc(6)

                // Inner loop

                Tic(7)

                for (int color = 0; color < 2; color++)
                {
                    forAllCells(xv, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xv(i,j,k) = bv(i,j,k);

                    forAllCells(xv, i, j, k)
                        if ((i + j + k) % 2 == color)
                            for (int c = 1; c < A.size(); c++)
                                xv(i,j,k) -=
                                    S(i,j,k)[c]
                                  * xv(labelVector(i,j,k) + offsets[c]);

                    forAllCells(xv, i, j, k)
                        if ((i + j + k) % 2 == color)
                            xv(i,j,k) /= S(i,j,k)[0];
                }

                Toc(7)

        // LEXGS (can only be done with an inner coefficient loop)

        // Scalar

            // SoA storage

                // Inner loop

                Tic(8)

                forAllCells(xs, i, j, k)
                {
                    xs(i,j,k) = bs(i,j,k);

                    for (int c = 1; c < A.size(); c++)
                        xs(i,j,k) -=
                            A[c](i,j,k)
                          * xs(labelVector(i,j,k) + offsets[c]);

                    xs(i,j,k) /= A[0](i,j,k);
                }

                Toc(8)

            // AoS storage

                // Inner loop

                Tic(9)

                forAllCells(xs, i, j, k)
                {
                    xs(i,j,k) = bs(i,j,k);

                    for (int c = 1; c < A.size(); c++)
                        xs(i,j,k) -=
                            S(i,j,k)[c]
                          * xs(labelVector(i,j,k) + offsets[c]);

                    xs(i,j,k) /= S(i,j,k)[0];
                }

                Toc(9)



        // Vector

            // SoA storage

                // Inner loop

                Tic(10)

                forAllCells(xs, i, j, k)
                {
                    xv(i,j,k) = bv(i,j,k);

                    for (int c = 1; c < A.size(); c++)
                        xv(i,j,k) -=
                            A[c](i,j,k)
                          * xv(labelVector(i,j,k) + offsets[c]);

                    xv(i,j,k) /= A[0](i,j,k);
                }

                Toc(10)

            // AoS storage

                // Inner loop

                Tic(11)

                forAllCells(xv, i, j, k)
                {
                    xv(i,j,k) = bv(i,j,k);

                    for (int c = 1; c < A.size(); c++)
                        xv(i,j,k) -=
                            S(i,j,k)[c]
                          * xv(labelVector(i,j,k) + offsets[c]);

                    xv(i,j,k) /= S(i,j,k)[0];
                }

                Toc(11)

        // JAC

        colocatedScalarField ys(xs);
        colocatedVectorField yv(xv);

        // Scalar

            // SoA storage

                // Outer loop

                Tic(12)

                forAllCells(xs, i, j, k)
                    ys(i,j,k) = bs(i,j,k);

                for (int c = 1; c < A.size(); c++)
                    forAllCells(xs, i, j, k)
                        ys(i,j,k) -=
                            A[c](i,j,k)
                          * xs(labelVector(i,j,k) + offsets[c]);

                forAllCells(xs, i, j, k)
                    ys(i,j,k) /= A[0](i,j,k);

                Toc(12)

                // Inner loop

                Tic(13)

                forAllCells(xs, i, j, k)
                    ys(i,j,k) = bs(i,j,k);

                forAllCells(xs, i, j, k)
                    for (int c = 1; c < A.size(); c++)
                        ys(i,j,k) -=
                            A[c](i,j,k)
                          * xs(labelVector(i,j,k) + offsets[c]);

                forAllCells(xs, i, j, k)
                    ys(i,j,k) /= A[0](i,j,k);

                Toc(13)

            // AoS storage

                // Outer loop

                Tic(14)

                forAllCells(xs, i, j, k)
                    ys(i,j,k) = bs(i,j,k);

                for (int c = 1; c < A.size(); c++)
                    forAllCells(S, i, j, k)
                        ys(i,j,k) -=
                            S(i,j,k)[c]
                          * xs(labelVector(i,j,k) + offsets[c]);

                forAllCells(xs, i, j, k)
                    ys(i,j,k) /= S(i,j,k)[0];

                Toc(14)

                // Inner loop

                Tic(15)

                forAllCells(xs, i, j, k)
                    ys(i,j,k) = bs(i,j,k);

                forAllCells(xs, i, j, k)
                    for (int c = 1; c < A.size(); c++)
                        ys(i,j,k) -=
                            S(i,j,k)[c]
                          * xs(labelVector(i,j,k) + offsets[c]);

                forAllCells(xs, i, j, k)
                    ys(i,j,k) /= S(i,j,k)[0];

                Toc(15)



        // Vector

            // SoA storage

                // Outer loop

                Tic(16)

                forAllCells(xv, i, j, k)
                    yv(i,j,k) = bv(i,j,k);

                for (int c = 1; c < A.size(); c++)
                    forAllCells(A[c], i, j, k)
                        yv(i,j,k) -=
                            A[c](i,j,k)
                          * xv(labelVector(i,j,k) + offsets[c]);

                forAllCells(xv, i, j, k)
                    yv(i,j,k) /= A[0](i,j,k);

                Toc(16)

                // Inner loop

                Tic(17)

                forAllCells(xv, i, j, k)
                    yv(i,j,k) = bv(i,j,k);

                forAllCells(xv, i, j, k)
                    for (int c = 1; c < A.size(); c++)
                        yv(i,j,k) -=
                            A[c](i,j,k)
                          * xv(labelVector(i,j,k) + offsets[c]);

                forAllCells(xv, i, j, k)
                    yv(i,j,k) /= A[0](i,j,k);

                Toc(17)

            // AoS storage

                // Outer loop

                Tic(18)

                forAllCells(xv, i, j, k)
                    yv(i,j,k) = bv(i,j,k);

                for (int c = 1; c < A.size(); c++)
                    forAllCells(S, i, j, k)
                        yv(i,j,k) -=
                            S(i,j,k)[c]
                          * xv(labelVector(i,j,k) + offsets[c]);

                forAllCells(xv, i, j, k)
                    yv(i,j,k) /= S(i,j,k)[0];

                Toc(18)

                // Inner loop

                Tic(19)

                forAllCells(xv, i, j, k)
                    yv(i,j,k) = bv(i,j,k);

                forAllCells(xv, i, j, k)
                    for (int c = 1; c < A.size(); c++)
                        yv(i,j,k) -=
                            S(i,j,k)[c]
                          * xv(labelVector(i,j,k) + offsets[c]);

                forAllCells(xv, i, j, k)
                    yv(i,j,k) /= S(i,j,k)[0];

                Toc(19)


        // Residual computation

        colocatedScalarField rs(xs);
        colocatedVectorField rv(xv);

        // Scalar

            // SoA storage

                // Outer loop

                Tic(20)

                forAllCells(xs, i, j, k)
                    rs(i,j,k) = bs(i,j,k);

                for (int c = 0; c < A.size(); c++)
                    forAllCells(xs, i, j, k)
                        rs(i,j,k) -=
                            A[c](i,j,k)
                          * xs(labelVector(i,j,k) + offsets[c]);

                Toc(20)

                // Inner loop

                Tic(21)

                forAllCells(xs, i, j, k)
                    rs(i,j,k) = bs(i,j,k);

                forAllCells(xs, i, j, k)
                    for (int c = 0; c < A.size(); c++)
                        rs(i,j,k) -=
                            A[c](i,j,k)
                          * xs(labelVector(i,j,k) + offsets[c]);

                Toc(21)

            // AoS storage

                // Outer loop

                Tic(22)

                forAllCells(xs, i, j, k)
                    rs(i,j,k) = bs(i,j,k);

                for (int c = 0; c < A.size(); c++)
                    forAllCells(S, i, j, k)
                        rs(i,j,k) -=
                            S(i,j,k)[c]
                          * xs(labelVector(i,j,k) + offsets[c]);

                Toc(22)

                // Inner loop

                Tic(23)

                forAllCells(xs, i, j, k)
                    rs(i,j,k) = bs(i,j,k);

                forAllCells(xs, i, j, k)
                    for (int c = 0; c < A.size(); c++)
                        rs(i,j,k) -=
                            S(i,j,k)[c]
                          * xs(labelVector(i,j,k) + offsets[c]);

                Toc(23)



        // Vector

            // SoA storage

                // Outer loop

                Tic(24)

                forAllCells(xv, i, j, k)
                    rv(i,j,k) = bv(i,j,k);

                for (int c = 0; c < A.size(); c++)
                    forAllCells(A[c], i, j, k)
                        rv(i,j,k) -=
                            A[c](i,j,k)
                          * xv(labelVector(i,j,k) + offsets[c]);

                Toc(24)

                // Inner loop

                Tic(25)

                forAllCells(xv, i, j, k)
                    rv(i,j,k) = bv(i,j,k);

                forAllCells(xv, i, j, k)
                    for (int c = 0; c < A.size(); c++)
                        rv(i,j,k) -=
                            A[c](i,j,k)
                          * xv(labelVector(i,j,k) + offsets[c]);

                Toc(25)

            // AoS storage

                // Outer loop

                Tic(26)

                forAllCells(xv, i, j, k)
                    rv(i,j,k) = bv(i,j,k);

                for (int c = 0; c < A.size(); c++)
                    forAllCells(S, i, j, k)
                        rv(i,j,k) -=
                            S(i,j,k)[c]
                          * xv(labelVector(i,j,k) + offsets[c]);

                Toc(26)

                // Inner loop

                Tic(27)

                forAllCells(xv, i, j, k)
                    rv(i,j,k) = bv(i,j,k);

                forAllCells(xv, i, j, k)
                    for (int c = 0; c < A.size(); c++)
                        rv(i,j,k) -=
                            S(i,j,k)[c]
                          * xv(labelVector(i,j,k) + offsets[c]);

                Toc(27)

        // Transfer from SoA to AoS

            // Outer loop

            Tic(28)
            colocatedStencilField S2("S2", fvMsh);
            forAll(A, c)
                forAllCells(A[c], i, j, k)
                    S2(i,j,k)[c] = A[c](i,j,k);
            Toc(28)

            // Inner loop

            Tic(29)
            colocatedStencilField S3("S3", fvMsh);
            forAllCells(S2, i, j, k)
                forAll(A, c)
                    S3(i,j,k)[c] = A[c](i,j,k);
            Toc(29)

        // Transfer from AoS to SoA

            // Outer loop

            Tic(30)
            colocatedScalarStencilField A2("A2", fvMsh);
            forAll(A, c)
                forAllCells(A[c], i, j, k)
                    A2[c](i,j,k) = S(i,j,k)[c];
            Toc(30)

            // Inner loop

            Tic(31)
            colocatedScalarStencilField A3("A3", fvMsh);
            forAllCells(S2, i, j, k)
                forAll(A, c)
                    A3[c](i,j,k) = S(i,j,k)[c];
            Toc(31)
    }

    Info<< "RBGS" << endl;

    int c = 0;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                Info<< (i == 0 ? "scalar " : "vector ")
                    << (j == 0 ? "SoA " : "AoS ")
                    << (k == 0 ? "outer " : "inner ")
                    << "= " << ticTocs[c++]/1000.0 << " ms"
                    << endl;

    Info<< endl;

    Info<< "LEXGS" << endl;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            Info<< (i == 0 ? "scalar " : "vector ")
                << (j == 0 ? "SoA " : "AoS ")
                << "= " << ticTocs[c++]/1000.0 << " ms"
                << endl;

    Info<< endl;

    Info<< "JAC" << endl;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                Info<< (i == 0 ? "scalar " : "vector ")
                    << (j == 0 ? "SoA " : "AoS ")
                    << (k == 0 ? "outer " : "inner ")
                    << "= " << ticTocs[c++]/1000.0 << " ms"
                    << endl;

    Info<< endl;

    Info<< "Residual" << endl;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                Info<< (i == 0 ? "scalar " : "vector ")
                    << (j == 0 ? "SoA " : "AoS ")
                    << (k == 0 ? "outer " : "inner ")
                    << "= " << ticTocs[c++]/1000.0 << " ms"
                    << endl;

    Info<< endl;

    Info<< "Storage format transfer" << endl;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            Info<< (i == 0 ? "SoA to AoS " : "AoS to SoA ")
                << (j == 0 ? "outer " : "inner ")
                << "= " << ticTocs[c++]/1000.0 << " ms"
                << endl;

    Info<< endl;
}
