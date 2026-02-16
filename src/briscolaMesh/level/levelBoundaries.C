#include "levelBoundaries.H"

#include "level.H"
#include "mesh.H"

#include "decompositionInterface.H"

#include "periodicPatch.H"
#include "patchBoundary.H"
#include "emptyBoundary.H"
#include "parallelBoundary.H"
#include "periodicBoundary.H"

namespace Foam
{

namespace briscola
{

void levelBoundaries::add
(
    const word name,
    const word type,
    const labelTensor T,
    const labelVector offset,
    const vector mappingOffset,
    const labelVector neighborOffset,
    const label neighborProcNum
)
{
    dictionary dict;

    dict.add("name", name);
    dict.add("type", type);
    dict.add("T", T);
    dict.add("offset", offset);

    if (mappingOffset != vector::zero)
        dict.add("mappingOffset", mappingOffset);

    if (neighborOffset != zeroXYZ)
        dict.add("neighborOffset", neighborOffset);

    if (neighborProcNum >= 0)
        dict.add("neighborProcNum", neighborProcNum);

    this->append(boundary::New(lvl_, dict));
}

void levelBoundaries::generateEmptyBoundaries()
{
    for (int f = 0; f < 6; f++)
        add
        (
            word
            (
                "empty-"
               + Foam::name(Pstream::myProcNo())
               + "-" + Foam::name(f)
            ),
            emptyBoundary::typeName,
            eye,
            faceOffsets[f]
        );

    // Also add dummy patch boundaries so that collective parallel operations
    // work. Not needed for periodic patches.

    forAll(lvl_.msh().patches(), patchi)
    if (!lvl_.msh().patches()[patchi].castable<periodicPatch>())
    {
        const patch& p = lvl_.msh().patches()[patchi];

        add
        (
            p.name(),
            p.type(),
            eye,
            zeroXYZ,
            p.dict().lookupOrDefault("offset",vector::zero)
        );
    }
}

void levelBoundaries::generateBrickInternalBoundaries()
{
    // Add parallel brick-internal boundaries

    const label myBrickNum = lvl_.decomp().myBrickNum();
    const labelVector myBrickDecomp = lvl_.decomp().myBrickDecomp();
    const labelVector myBrickPart = lvl_.decomp().myBrickPart();

    if (cmptProduct(myBrickDecomp) > 1)
    {
        // Try all boundary offset vectors. Even if the bricks are unstructured,
        // within a brick we can add edge and vertex boundaries anyway.

        const labelVector S =
            cmptMax(myBrickPart - unitXYZ, zeroXYZ);

        const labelVector E =
            cmptMin(myBrickPart + unitXYZ, myBrickDecomp - unitXYZ);

        labelVector nei;
        for (nei.x() = S.x(); nei.x() <= E.x(); nei.x()++)
        for (nei.y() = S.y(); nei.y() <= E.y(); nei.y()++)
        for (nei.z() = S.z(); nei.z() <= E.z(); nei.z()++)
        if (nei != myBrickPart)
        {
            const labelVector bo = nei - myBrickPart;

            const label neighborProcNum =
                lvl_.decomp().brickProcMaps()[myBrickNum](nei);

            // Add to boundaries

            add
            (
                word
                (
                    "parallel-" + Foam::name(Pstream::myProcNo())
                  + "to" + Foam::name(neighborProcNum)
                ),
                "parallel",
                eye,
                bo,
                vector::zero,
              - bo,
                neighborProcNum
            );
        }
    }
}

void levelBoundaries::generateBrickExternalBoundaries()
{
    const brickTopology& topo = lvl_.msh().topology();

    forAll(lvl_.msh().bricks(), bricki)
    {
        forAll(topo.links()[bricki].faceLinks(), linki)
        if (topo.links()[bricki].faceLinks().set(linki))
        {
            const brickLink link(topo.links()[bricki].faceLinks()[linki]);

            generateLinkBoundaries(link);
        }

        forAll(topo.links()[bricki].edgeLinks(), linki)
        if (topo.links()[bricki].edgeLinks().set(linki))
        {
            const brickLink link(topo.links()[bricki].edgeLinks()[linki]);

            generateLinkBoundaries(link);
        }

        forAll(topo.links()[bricki].vertexLinks(), linki)
        if (topo.links()[bricki].vertexLinks().set(linki))
        {
            const brickLink link(topo.links()[bricki].vertexLinks()[linki]);

            generateLinkBoundaries(link);
        }
    }
}

void levelBoundaries::generateLinkBoundaries(const brickLink& link)
{
    const decompositionInterface interface(link, lvl_);

    const labelBlock& map = interface.map();
    const PtrList<labelPair>& slices = interface.slices();

    const labelVector bo = link.offset();
    const labelTensor T = link.T();

    // Only add slices if they belong to this processor

    forAllBlock(map, i, j, k)
    if (slices[map(i,j,k)].first() == Pstream::myProcNo())
    {
        const labelVector ijk(i,j,k);
        const labelPair& slice = slices[map(ijk)];

        add
        (
            word
            (
                word(link.periodic() ? "periodic" : "parallel")
              + "-" + Foam::name(slice.first())
              + "to" + Foam::name(slice.second())
            ),
            link.periodic() ? "periodic" : "parallel",
            T,
            bo,
            vector::zero,
          - (T.T() & bo),
            slice.second()
        );

        // Add internal edges and vertices of face boundaries

        if (cmptSum(cmptMag(bo)) == 1)
        {
            const label faceNormalDir = faceNumber(bo)/2;

            // Test in the two directions perpendicular to the face normal

            const labelVector dir1 = (faceNormalDir == 0 ? unitY : unitX);
            const labelVector dir2 = (faceNormalDir == 2 ? unitY : unitZ);

            for (label ii = -1; ii <= 1; ii++)
            for (label jj = -1; jj <= 1; jj++)
            if (ii != 0 || jj != 0)
            {
                const labelVector ijk2(ijk + dir1*ii + dir2*jj);
                const labelVector bo2(bo + dir1*ii + dir2*jj);

                if
                (
                    cmptMin(ijk2) >= 0
                 && ijk2.x() < map.l()
                 && ijk2.y() < map.m()
                 && ijk2.z() < map.n()
                )
                {
                    const labelPair& slice2 = slices[map(ijk2)];

                    add
                    (
                        word
                        (
                            word(link.periodic() ? "periodic" : "parallel")
                          + "-" + Foam::name(slice.first())
                          + "to" + Foam::name(slice2.second())
                        ),
                        link.periodic() ? "periodic" : "parallel",
                        T,
                        bo2,
                        vector::zero,
                      - (T.T() & bo2),
                        slice2.second()
                    );
                }
            }
        }

        // Add internal vertices of edge boundaries

        if (cmptSum(cmptMag(bo)) == 2)
        {
            const label edgeDir = edgeNumber(bo)/4;

            // Test in the direction of the edge

            const labelVector dir = units[edgeDir];

            for (label ii = -1; ii <= 1; ii+=2)
            {
                const labelVector ijk2(ijk + dir*ii);
                const labelVector bo2(bo + dir*ii);

                if
                (
                    cmptMin(ijk2) >= 0
                 && ijk2.x() < map.l()
                 && ijk2.y() < map.m()
                 && ijk2.z() < map.n()
                )
                {
                    const labelPair& slice2 = slices[map(ijk2)];

                    add
                    (
                        word
                        (
                            word(link.periodic() ? "periodic" : "parallel")
                          + "-" + Foam::name(slice.first())
                          + "to" + Foam::name(slice2.second())
                        ),
                        link.periodic() ? "periodic" : "parallel",
                        T,
                        bo2,
                        vector::zero,
                      - (T.T() & bo2),
                        slice2.second()
                    );
                }
            }
        }
    }
}

void levelBoundaries::generatePatchBoundaries()
{
    FastPtrList<boundary>& boundaries = *this;

    // Loop over all patches and find the corresponding brick face. If not
    // found, the faces of this level are not part of the patch and we set a
    // fake patch boundary associated with the patch, with a zero boundary
    // offset. This is needed to enable global parallel operations on patch
    // boundaries.

    const label brickNum = lvl_.decomp().myBrickNum();
    const brick& b = lvl_.msh().bricks()[brickNum];

    forAll(lvl_.msh().patches(), patchi)
    if (!lvl_.msh().patches()[patchi].castable<periodicPatch>())
    {
        const patch& p = lvl_.msh().patches()[patchi];

        bool found = false;

        // Multiple faces can belong to the same boundary, so don't break the
        // loop after finding a matching face.

        forAll(b.faces(), facei)
        {
            const face& fi = b.faces()[facei];

            // Skip faces that are already assigned to a boundary

            bool assigned = false;
            forAll(boundaries, i)
                if (boundaries[i].offset() == faceOffsets[facei])
                    assigned = true;

            if (assigned)
                continue;

            forAll(p.facePtrs(), facej)
            {
                const face& fj = *p.facePtrs()[facej];

                if (&fi == &fj)
                {
                    add
                    (
                        p.name(),
                        p.type(),
                        eye,
                        faceOffsets[facei],
                        p.dict().lookupOrDefault("offset", vector::zero)
                    );

                    found = true;
                }
            }
        }

        // Set fake patch with zero boundary offset

        if (!found)
        {
            add
            (
                p.name(),
                p.type(),
                eye,
                zeroXYZ,
                p.dict().lookupOrDefault("offset",vector::zero)
            );
        }
    }
}

void levelBoundaries::checkFaceBoundaries()
{
    const FastPtrList<boundary>& boundaries = *this;

    for (int f = 0; f < 6; f++)
    {
        label count = 0;
        forAll(boundaries, j)
            if (faceOffsets[f] == boundaries[j].offset())
                count++;

        if (count < 1)
            FatalErrorInFunction
                << "Could not find the patch to which face " << f
                << " belongs on level " << lvl_.levelNum() << endl
                << abort(FatalError);

        if (count > 1)
            FatalErrorInFunction
                << "Face " << f << " is associated with "
                << count << " boundaries on level " << lvl_.levelNum()
                << endl << abort(FatalError);
    }
}

void levelBoundaries::setCommTags()
{
    FastPtrList<boundary>& boundaries = *this;

    const label n = Pstream::nProcs();
    const label m = Pstream::myProcNo();

    List<List<dictionary>> bDicts(n, List<dictionary>(0));

    forAll(boundaries, bi)
        if (boundaries[bi].castable<parallelBoundary>())
            bDicts[m].append(dictionary(boundaries[bi].dict()));

    Pstream::gatherList(bDicts);

    labelList myTags;

    if (Pstream::master())
    {
        List<labelList> tags(n,labelList(0));

        forAll(tags, i)
        {
            tags[i].setSize(bDicts[i].size());
            tags[i] = 1;
        }

        // The pairCount array can get pretty big for a large number of
        // processors, but allocating an n by n array is by far the easiest way
        // to track pairs. On 10.000 processors, the array will consume 381 MB,
        // which should be doable.

        List<labelList> pairCount(n, labelList(n));

        for (label i = 0; i < n; i++)
            for (label j = 0; j < n; j++)
                pairCount[i][j] = 0;

        forAll(bDicts, i)
        forAll(bDicts[i], bi)
        {
            const dictionary& dict0 = bDicts[i][bi];

            const labelVector source
            (
              - labelTensor(dict0.lookup("T")).T()
              & labelVector(dict0.lookup("offset"))
            );

            const label j = readLabel(dict0.lookup("neighborProcNum"));

            // Avoid searching pairs twice

            if (i <= j)
            {
                bool found = false;

                forAll(bDicts[j], bj)
                if (i != j || bi != bj)
                {
                    const dictionary& dict1 = bDicts[j][bj];

                    const label k =
                        readLabel(dict1.lookup("neighborProcNum"));

                    const labelVector target =
                        dict1.lookup("offset");

                    if (k == i && source == target)
                    {
                        pairCount[i][j]++;

                        tags[i][bi] = pairCount[i][j];
                        tags[j][bj] = pairCount[i][j];

                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    FatalErrorInFunction
                        << "Could not find neighbor of boundary "
                        << word(dict0.lookup("name"))
                        << " of processor " << j << endl
                        << abort(FatalError);
                }
            }
        }

        // Send

        forAll(bDicts, i)
        {
            if (i == Pstream::masterNo())
            {
                myTags = tags[i];
            }
            else
            {
                OPstream toSlave
                (
                    Pstream::commsTypes::blocking,
                    i
                );

                toSlave << tags[i];
            }
        }
    }
    else
    {
        // Receive

        IPstream fromMaster
        (
            Pstream::commsTypes::blocking,
            Pstream::master()
        );

        fromMaster >> myTags;
    }

    label c = 0;

    forAll(boundaries, bi)
    {
        boundary& b = boundaries[bi];

        if (b.castable<parallelBoundary>())
            b.cast<parallelBoundary>().setTag(myTags[c++]);
    }
}

void levelBoundaries::setMask()
{
    FastPtrList<boundary>& boundaries = *this;

    mask_.setSize(3,3,3);
    parallelMask_.setSize(3,3,3);

    mask_ = Zero;
    parallelMask_ = Zero;

    // Set boundary mask

    forAll(boundaries, i)
        if (boundaries[i].offsetDegree() > 0)
            mask_(boundaries[i].offset()+unitXYZ) = 1;

    // Set parallel/periodic boundary mask

    forAll(boundaries, i)
        if (boundaries[i].offsetDegree() > 0)
            if (boundaries[i].castable<parallelBoundary>())
                parallelMask_(boundaries[i].offset()+unitXYZ) = 1;
}

void levelBoundaries::reorder()
{
    FastPtrList<boundary>& boundaries = *this;

    labelList oldToNew(boundaries.size());

    int j = 0;

    // Patch boundaries

    forAll(boundaries, i)
        if (boundaries[i].castable<patchBoundary>())
            oldToNew[i] = j++;

    // Parallel/periodic boundaries

    for (int d = 1; d <= 3; d++)
        forAll(boundaries, i)
            if
            (
                boundaries[i].castable<parallelBoundary>()
             && boundaries[i].offsetDegree() == d
            )
                oldToNew[i] = j++;

    // Empty boundaries

    forAll(boundaries, i)
        if (boundaries[i].castable<emptyBoundary>())
            oldToNew[i] = j++;

    if (j != boundaries.size())
        FatalErrorInFunction
            << "Could not order boundaries"
            << endl << abort(FatalError);

    boundaries.reorder(oldToNew);
}

levelBoundaries::levelBoundaries(const level& lvl)
:
    FastPtrList<boundary>(0),
    lvl_(lvl)
{
    if (lvl_.empty())
    {
        generateEmptyBoundaries();
    }
    else
    {
        generateBrickInternalBoundaries();
        generateBrickExternalBoundaries();
        generatePatchBoundaries();
        reorder();
    }

    checkFaceBoundaries();

    setCommTags();
    setMask();
}

levelBoundaries::levelBoundaries(const levelBoundaries& b)
:
    FastPtrList<boundary>(b),
    lvl_(b.lvl_),
    mask_(b.mask_),
    parallelMask_(b.parallelMask_)
{}

levelBoundaries::levelBoundaries(const levelBoundaries& b, const level& lvl)
:
    FastPtrList<boundary>(b, lvl),
    lvl_(lvl),
    mask_(b.mask_),
    parallelMask_(b.parallelMask_)
{}

levelBoundaries::~levelBoundaries()
{}

}

}
