#include "decompositionMap.H"
#include "decomposition.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

decompositionMap::decompositionMap(const decomposition& decomp)
:
    decomp_(decomp)
{
    // Determine map shape

    const labelBlock& brickMap = decomp_.lvl().msh().topology().map();

    PtrList<labelList> N(3);
    PtrList<labelList> M(3);

    forAll(N, dir)
    {
        N.set(dir, new labelList(brickMap.shape()[dir], 0));
        M.set(dir, new labelList(brickMap.shape()[dir], 0));
    }

    forAllBlock(brickMap, i, j, k)
    {
        const labelVector ijk(i,j,k);
        const label brick = brickMap(ijk);

        if (brick > -1)
        {
            const labelVector D = decomp_.brickDecomps()[brick];
            const labelVector P =
                cmptDivide(decomp_.lvl().msh().bricks()[brick].N(), D);

            for (int dir = 0; dir < 3; dir++)
            {
                const label l = ijk[dir];

                if (N[dir][l] == 0)
                {
                    N[dir][l] = D[dir];
                    M[dir][l] = P[dir];
                }
            }
        }
    }

    const labelVector S(sum(N[0]), sum(N[1]), sum(N[2]));

    // Set map

    labelBlock& map = *this;

    map.setSize(S);
    map = -1;

    labelVector cursor(zeroXYZ);

    for (int i = 0; i < brickMap.l(); i++)
    {
        cursor.y() = 0;

        for (int j = 0; j < brickMap.m(); j++)
        {
            cursor.z() = 0;

            for (int k = 0; k < brickMap.n(); k++)
            {
                const labelVector ijk(i,j,k);
                const label brick = brickMap(ijk);

                // Set processor numbers

                if (brick > -1)
                    forAllBlock(decomp_.brickProcMaps()[brick], a, b, c)
                        map(cursor + labelVector(a,b,c)) =
                            decomp_.brickProcMaps()[brick](a,b,c);

                cursor.z() += N[2][k];
            }

            cursor.y() += N[1][j];
        }

        cursor.x() += N[0][i];
    }

    // Set legend

    legend_.clear();

    forAllBlock(map, i, j, k)
    if (map(i,j,k) != -1)
    {
        label v = map(i,j,k);

        if (v >= legend_.size())
            legend_.setSize(v+1);

        legend_[v] = labelVector(i,j,k);
    }

    // Set shapes and starts

    shapes_.setSize(S);
    starts_.setSize(S);
    starts_ = zeroXYZ;

    cursor = zeroXYZ;

    for (int i = 0; i < brickMap.l(); i++)
    {
        cursor.y() = 0;

        for (int j = 0; j < brickMap.m(); j++)
        {
            cursor.z() = 0;

            for (int k = 0; k < brickMap.n(); k++)
            {
                const labelVector ijk(i,j,k);

                // Part shape in this brick

                const labelVector shape(M[0][i], M[1][j], M[2][k]);

                // Brick decomposition shape

                const labelVector brickDecompShape(N[0][i], N[1][j], N[2][k]);

                // Start index of current brick

                labelVector start;

                for (int dir = 0; dir < 3; dir++)
                    start[dir] =
                        ijk[dir] > 0
                      ? starts_
                        (
                            legend_
                            (
                                decomp_.brickProcMaps()
                                [
                                    brickMap(ijk-units[dir])
                                ](zeroXYZ)
                            )
                        )[dir]
                      + decomp_.brickDecomps()[brickMap(ijk-units[dir])][dir]
                      * shapes_
                        (
                            legend_
                            (
                                decomp_.brickProcMaps()
                                [
                                    brickMap(ijk-units[dir])
                                ](zeroXYZ)
                            )
                        )[dir]
                      : 0;

                // Set shapes and starts

                for (int a = 0; a < brickDecompShape.x(); a++)
                for (int b = 0; b < brickDecompShape.y(); b++)
                for (int c = 0; c < brickDecompShape.z(); c++)
                {
                    const labelVector abc(a,b,c);

                    shapes_(cursor + abc) = shape;
                    starts_(cursor + abc) = start + cmptMultiply(abc, shape);
                }

                cursor.z() += N[2][k];
            }

            cursor.y() += N[1][j];
        }

        cursor.x() += N[0][i];
    }
}

decompositionMap::decompositionMap(const decompositionMap& map)
:
    labelBlock(map),
    decomp_(map.decomp_),
    legend_(map.legend_),
    shapes_(map.shapes_),
    starts_(map.starts_)
{}

decompositionMap::~decompositionMap()
{}

}

}
