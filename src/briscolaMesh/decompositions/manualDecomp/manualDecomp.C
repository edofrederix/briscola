#include "manualDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

namespace decompositions
{

defineTypeNameAndDebug(manualDecomp, 0);
addToRunTimeSelectionTable(decomposition, manualDecomp, dictionary);

manualDecomp::manualDecomp(mesh& msh)
:
    decomposition(msh)
{
    List<labelVector> brickDecomps(dict_.lookup("brickDecompositions"));

    // Check if all bricks are specified

    if (brickDecomps.size() != msh.bricks().size())
        FatalErrorInFunction
            << "The number of manual brick decomposition vectors should be "
            << "equal to the number of bricks" << abort(FatalError);

    // Check if the total number of processors matches

    label nProcs(0);

    forAll(brickDecomps, bricki)
    {
        const label nProcsi = cmptProduct(brickDecomps[bricki]);

        if (nProcsi < 1)
            FatalErrorInFunction
                << "Incorrect decomposition given for " << bricki << endl
                << abort(FatalError);

        nProcs += nProcsi;
    }

    if (nProcs != Pstream::nProcs())
        FatalErrorInFunction
            << "Mismatch between the specified number of processors ("
            << nProcs << ") and the actual number of processors ("
            << Pstream::nProcs() << ")" << endl;

    // Transform decompositions according to the brick's transformation

    forAll(brickDecomps, i)
        brickDecomps[i] = cmptMag(msh.bricks()[i].T() & brickDecomps[i]);

    // Check if the decompositions are feasible

    forAll(msh.bricks(), i)
    {
        const labelVector N = msh.bricks()[i].N();
        const labelVector D = brickDecomps[i];

        for (label dir = 0; dir < 3; dir++)
        {
            if (N[dir] % D[dir] != 0)
                FatalErrorInFunction
                    << "Brick " << i << " has " << N[dir] << " cells in the "
                    << dir << " direction but this is incompatible with a "
                    << "decomposition of " << D[dir] << " processors along "
                    << "this direction" << endl << abort(FatalError);

            label Nd = N[dir]/D[dir];

            if (!Nd || ((Nd & (Nd-1)) != 0 && ((Nd/3) & ((Nd/3)-1)) != 0))
                FatalErrorInFunction
                    << "The number of cells (" << Nd << ") of the mesh part "
                    << "that results from the decomposition of brick " << i
                    << " in the " << dir << " direction is not a power of 2 "
                    << "nor a triple of a power of 2" << endl
                    << abort(FatalError);
        }
    }

    // Check if decompositions agree across brick face links

    const brickTopology& topo = msh.topology();

    forAll(msh.bricks(), bricki)
    {
        const brickLinks& links = topo.links()[bricki];

        forAll(links.faceLinks(), i)
        if (links.faceLinks().set(i))
        {
            const brickFaceLink& link = links.faceLinks()[i];

            const label brickj = link.f1().parentBrick().num();

            const label facei = link.f0().num();
            const label facej = link.f1().num();

            labelVector Ni = brickDecomps[bricki];
            labelVector Nj = brickDecomps[brickj];

            Ni[facei/2] = 1;
            Nj[facej/2] = 1;

            Nj = cmptMag(link.T() & Nj);

            if (Ni != Nj && Pstream::master())
                FatalErrorInFunction
                    << "Inconsistent decomposition between brick " << bricki
                    << " on face " << facei << " and brick " << brickj
                    << " on face " << facej << endl
                    << abort(FatalError);
        }
    }

    // Distribute brick parts

    label proc(0);

    forAll(msh.bricks(), bricki)
    {
        const labelVector& decomp = brickDecomps[bricki];

        for (label i = 0; i < decomp.x(); i++)
        {
            for (label j = 0; j < decomp.y(); j++)
            {
                for (label k = 0; k < decomp.z(); k++)
                {
                    if (proc == Pstream::myProcNo())
                    {
                        myBrickNum_ = bricki;
                        myBrickDecomp_ = decomp;
                        myBrickPart_ = labelVector(i,j,k);
                    }

                    proc++;
                }
            }
        }
    }

    init();
}

manualDecomp::manualDecomp(const manualDecomp& decomp)
:
    decomposition(decomp),
    myBrickNum_(decomp.myBrickNum_),
    myBrickDecomp_(decomp.myBrickDecomp_),
    myBrickPart_(decomp.myBrickPart_)
{}

}

}

}