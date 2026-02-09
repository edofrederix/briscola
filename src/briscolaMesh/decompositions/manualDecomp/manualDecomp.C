#include "manualDecomp.H"
#include "addToRunTimeSelectionTable.H"
#include "level.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

namespace decompositions
{

defineTypeNameAndDebug(manualDecomp, 0);
addToRunTimeSelectionTable(decomposition, manualDecomp, dictionary);

manualDecomp::manualDecomp(const level& lvl, const bool coarsen)
:
    decomposition(lvl)
{
    List<labelVector> brickDecomps(dict_.lookup("brickDecompositions"));

    // Check if all bricks are specified

    if (brickDecomps.size() != lvl.msh().bricks().size())
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
            << Pstream::nProcs() << ")" << endl << abort(FatalError);

    // Check if decomposition components are factors of two or a triple factor
    // of two

    forAll(brickDecomps, bricki)
    {
        for (label d = 0; d < 3; d++)
        {
            label Nd = brickDecomps[bricki][d];

            if (!powerOfTwo(Nd) && ((Nd % 3) || !powerOfTwo(Nd/3)))
                FatalErrorInFunction
                    << "The decomposition of brick " << bricki << "in the "
                    << d << " direction is not a power of 2 nor a triple of "
                    << "a power of 2" << endl << abort(FatalError);
        }
    }

    // Transform decompositions according to the brick's transformation

    forAll(brickDecomps, i)
        brickDecomps[i] = cmptMag(lvl.msh().bricks()[i].T() & brickDecomps[i]);

    // Check if the decompositions are feasible

    forAll(lvl.msh().bricks(), i)
    {
        const labelVector N = lvl.msh().bricks()[i].N();
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

            if (!powerOfTwo(Nd) && ((Nd % 3) || !powerOfTwo(Nd/3)))
                FatalErrorInFunction
                    << "The number of cells (" << Nd << ") of the brick part "
                    << "that results from the decomposition of brick " << i
                    << " in the " << dir << " direction is not a power of 2 "
                    << "nor a triple of a power of 2" << endl
                    << abort(FatalError);
        }
    }

    // Check if decompositions agree across brick face links

    const brickTopology& topo = lvl.msh().topology();

    forAll(lvl.msh().bricks(), bricki)
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

    // Coarsen the brick decomposition if needed, which is the case when the
    // parent level is not coarsenable

    const List<labelVector> brickDecomps0(brickDecomps);

    if (lvl_.hasParent() && !lvl_.parent().coarsenable())
    {
        // Use parent brick decompositions

        brickDecomps = lvl_.parent().decomp().brickDecomps();

        // Coarsen brick decomposition components. If a value is 3, then the
        // division by 2 gives 1 through integer cast.

        forAll(brickDecomps, bricki)
            for (int d = 0; d < 3; d++)
                brickDecomps[bricki][d] =
                    Foam::max(brickDecomps[bricki][d]/2, 1);

        // Check coarsened decompositions: the coarsened brick decomposition
        // needs to be an integer fraction of the set brick decomposition

        forAll(brickDecomps, bricki)
            for (int d = 0; d < 3; d++)
                if (brickDecomps0[bricki][d] % brickDecomps[bricki][d])
                    FatalErrorInFunction
                        << "Invalid decomposition coarsening"
                        << abort(FatalError);
    }

    // Distribute brick parts

    myBrickPart_ = -unitXYZ;

    label proc = 0;
    forAll(lvl.msh().bricks(), bricki)
    {
        const labelVector& decomp = brickDecomps[bricki];
        const labelVector& decomp0 = brickDecomps0[bricki];
        const labelVector R = cmptDivide(decomp0, decomp);

        // Create map without coarsening

        labelBlock map(decomp0);
        labelVector ijk;

        for (ijk.x() = 0; ijk.x() < decomp0.x(); ijk.x()++)
        for (ijk.y() = 0; ijk.y() < decomp0.y(); ijk.y()++)
        for (ijk.z() = 0; ijk.z() < decomp0.z(); ijk.z()++)
        {
            map(ijk) = proc;

            if (proc == Pstream::myProcNo())
            {
                myBrickNum_ = bricki;
                myBrickDecomp_ = decomp;
            }

            proc++;
        }

        // Sample map, possibly with coarsening

        for (ijk.x() = 0; ijk.x() < decomp.x(); ijk.x()++)
        for (ijk.y() = 0; ijk.y() < decomp.y(); ijk.y()++)
        for (ijk.z() = 0; ijk.z() < decomp.z(); ijk.z()++)
        {
            const labelVector ijkR(cmptMultiply(ijk, R));

            if (map(ijkR) == Pstream::myProcNo())
                myBrickPart_ = ijk;
        }
    }
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