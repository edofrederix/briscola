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
    decompPerBrick_ =
        List<labelVector>(dict_.lookup("brickDecompositions"));

    // Check if all bricks are specified

    if (decompPerBrick_.size() != msh.bricks().size())
    {
        FatalError
            << "The number of manualDecomp vectors should be equal "
            << "to the number of bricks" << endl;
        FatalError.exit();
    }

    // Check if the total number of processors matches

    label nProcs(0);

    forAll(decompPerBrick_, bricki)
    {
        const label nProcsi = cmptProduct(decompPerBrick_[bricki]);

        if (nProcsi < 1)
        {
            FatalError
                << "Incorrect decomposition given for " << bricki << endl;
            FatalError.exit();
        }

        nProcs += nProcsi;
    }

    if (nProcs != Pstream::nProcs())
    {
        FatalError
            << "Mismatch between the specified number of processors ("
            << nProcs << ") and the actual number of processors ("
            << Pstream::nProcs() << ")" << endl;
        FatalError.exit();
    }

    // Transform decompositions according to the brick's transformation

    forAll(decompPerBrick_, i)
    {
        const brick& b = msh.bricks()[i];

        decompPerBrick_[i] = cmptMag(b.T() & decompPerBrick_[i]);
    }

    // Check if the decompositions are feasible

    forAll(msh.bricks(), bricki)
    {
        const brick& b = msh.bricks()[bricki];

        const labelVector N = b.N();
        const labelVector D = decompPerBrick_[bricki];

        for (label dir = 0; dir < 3; dir++)
        {
            if (N[dir] % D[dir] != 0)
            {
                FatalError
                    << "Brick " << bricki << " has " << N[dir]
                    << " cells in the " << dir
                    << " direction but this is incompatible with a decomposition of "
                    << D[dir] << " processors along this direction"
                    << endl;
                FatalError.exit();
            }

            label Nd = N[dir]/D[dir];

            if (!Nd || ((Nd & (Nd-1)) != 0 && ((Nd/3) & ((Nd/3)-1)) != 0))
            {
                FatalError
                    << "The number of cells (" << Nd << ") of "
                    << "the mesh part that results from the decomposition of brick "
                    << bricki << " in the " << dir << " direction "
                    << "is not a power of 2 nor a triple of a power of 2" << endl;
                FatalError.exit();
            }
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

            labelVector Ni = decompPerBrick_[bricki];
            labelVector Nj = decompPerBrick_[brickj];

            Ni[facei/2] = 1;
            Nj[facej/2] = 1;

            Nj = cmptMag(link.T() & Nj);

            if (Ni != Nj && Pstream::master())
            {
                FatalErrorInFunction
                    << "Inconsistent decomposition between brick " << bricki
                    << " on face " << facei << " and brick " << brickj
                    << " on face " << facej << endl
                    << abort(FatalError);
            }
        }
    }

    // Distribute brick parts

    label proc(0);

    forAll(msh.bricks(), bricki)
    {
        const labelVector& decomp = decompPerBrick_[bricki];

        for (label i = 0; i < decomp.x(); i++)
        {
            for (label j = 0; j < decomp.y(); j++)
            {
                for (label k = 0; k < decomp.z(); k++)
                {
                    if (proc == Pstream::myProcNo())
                    {
                        myBrickNum_ = bricki;
                        myBrickPart_ = labelVector(i,j,k);
                    }

                    proc++;
                }
            }
        }
    }

    updateGlobalData(msh);
}

manualDecomp::manualDecomp(const manualDecomp& decomp)
:
    decomposition(decomp),
    decompPerBrick_(decomp.decompPerBrick_),
    myBrickNum_(decomp.myBrickNum_),
    myBrickPart_(decomp.myBrickPart_)
{}

}

}

}