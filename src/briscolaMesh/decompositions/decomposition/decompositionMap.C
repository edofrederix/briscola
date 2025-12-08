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

    const labelBlock& brickMap = decomp_.msh().topology().map();

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
                cmptDivide(decomp_.msh().bricks()[brick].N(), D);

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

    // Set map, shapes and starts

    labelBlock& map = *this;

    map.setSize(S);
    map = -1;

    shapes_.setSize(S);
    starts_.setSize(S);
    starts_ = zeroXYZ;

    labelVector cursor(zeroXYZ);

    for (int i = 0; i < brickMap.l(); i++)
    {
        cursor.y() = 0;

        for (int j = 0; j < brickMap.m(); j++)
        {
            cursor.z() = 0;

            for (int k = 0; k < brickMap.n(); k++)
            {
                const label brick = brickMap(i,j,k);

                // Set processor numbers

                if (brick > -1)
                    forAllBlock(decomp_.brickProcMaps()[brick], a, b, c)
                        map(cursor + labelVector(a,b,c)) =
                            decomp_.brickProcMaps()[brick](a,b,c);

                // Set shapes and starts

                for (int a = 0; a < N[0][i]; a++)
                for (int b = 0; b < N[1][j]; b++)
                for (int c = 0; c < N[2][k]; c++)
                {
                    const labelVector ijk(cursor + labelVector(a,b,c));

                    shapes_(ijk) = labelVector(M[0][i], M[1][j], M[2][k]);

                    for (int dir = 0; dir < 3; dir++)
                        if (ijk[dir]-1 >= 0)
                            starts_(ijk)[dir] =
                                starts_(ijk[dir]-1)[dir]
                              + shapes_(ijk[dir]-1)[dir];
                }

                cursor.z() += N[2][k];
            }

            cursor.y() += N[1][j];
        }

        cursor.x() += N[0][i];
    }

    // Legend

    legend_.clear();

    forAllBlock(map, i, j, k)
    if (map(i,j,k) != -1)
    {
        label v = map(i,j,k);

        if (v >= legend_.size())
            legend_.setSize(v+1);

        legend_[v] = labelVector(i,j,k);
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
