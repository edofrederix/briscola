#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(mesh, 0);

void mesh::generateBrickInternalPartPatches()
{
    // Add parallel brick-internal patches

    const labelVector myBrickDecomp(decomp().myBrickDecomp());
    const labelVector myBrickPart(decomp().myBrickPart());

    if (cmptProduct(myBrickDecomp) > 1)
    {
        // Try all boundary offset vectors, including faces, edges and vertices.
        // Even if the bricks are unstructured, within a brick we can add edge
        // and vertex part patches anyway.

        labelVector bo;

        for (bo.x() = -1; bo.x() < 2; bo.x()++)
        for (bo.y() = -1; bo.y() < 2; bo.y()++)
        for (bo.z() = -1; bo.z() < 2; bo.z()++)
        if (cmptSum(cmptMag(bo)) != 0)
        {
            labelVector neighbor(myBrickPart + bo);

            const label neighborProcNum
            (
                decomp().getProcNum
                (
                    neighbor,
                    decomp().myBrickNum()
                )
            );

            if (neighborProcNum == -1) continue;

            // Add to part patches

            dictionary dict;

            dict.add
            (
                "name",
                "parallel-"
              + Foam::name(Pstream::myProcNo())
              + "to"
              + Foam::name(neighborProcNum)
            );

            dict.add("type", "parallel");
            dict.add("T", eye);
            dict.add("boundaryOffset", bo);
            dict.add("neighborOffset", -bo);
            dict.add("neighborProcNum", neighborProcNum);

            partPatches_.append
            (
                partPatch::New(*this, dict)
            );
        }
    }
}

void mesh::generateBrickExternalPartPatches()
{
    const brickTopology& topo = topology();

    forAll(bricks(), bricki)
    {
        forAll(topo.links()[bricki].faceLinks(), linki)
        if (topo.links()[bricki].faceLinks().set(linki))
        {
            const brickLink link(topo.links()[bricki].faceLinks()[linki]);

            generateLinkPartPatches(link);
        }

        forAll(topo.links()[bricki].edgeLinks(), linki)
        if (topo.links()[bricki].edgeLinks().set(linki))
        {
            const brickLink link(topo.links()[bricki].edgeLinks()[linki]);

            generateLinkPartPatches(link);
        }

        forAll(topo.links()[bricki].vertexLinks(), linki)
        if (topo.links()[bricki].vertexLinks().set(linki))
        {
            const brickLink link(topo.links()[bricki].vertexLinks()[linki]);

            generateLinkPartPatches(link);
        }
    }
}

void mesh::generateLinkPartPatches(const brickLink& link)
{
    const brickDecompositionInterface interface(link, *this);

    const PtrList<brickDecompositionSlice>& slices = interface.slices();

    const labelVector bo = link.offset();
    const labelTensor T = link.T();

    // Only add slices if they belong to this processor

    forAll(slices, slicei)
    if (slices[slicei].procNum0() == Pstream::myProcNo())
    {
        const brickDecompositionSlice& slice = slices[slicei];

        dictionary dict;

        dict.add
        (
            "name",
            word(link.periodic() ? "periodic" : "parallel")
          + "-"
          + Foam::name(slice.procNum0())
          + "to"
          + Foam::name(slice.procNum1())
        );

        dict.add("type", link.periodic() ? "periodic" : "parallel");
        dict.add("T", T);
        dict.add("boundaryOffset", bo);
        dict.add("neighborOffset", -(T.T() & bo));
        dict.add("neighborProcNum", slice.procNum1());

        partPatches_.append(partPatch::New(*this, dict));
    }
}

void mesh::generateBoundaryPartPatches()
{
    // All face part patches that are not yet set must be part of a boundary
    // patch

    for (label facei = 0; facei < 6; facei++)
    {
        const labelVector offset = faceOffsets[facei];
        const label brickNum = decomp().myBrickNum();
        const brick& b = bricks()[brickNum];
        const face& f = b.faces()[facei];

        forAll(partPatches_, patchi)
        if (partPatches_[patchi].boundaryOffset() == offset)
        {
            goto found;
        }

        forAll(patches(), patchi)
        if (patches()[patchi].type() != patch::PERIODIC)
        {
            const patch& p = patches()[patchi];

            forAll(p.facePtrs(), facej)
            if (&f == p.facePtrs()[facej])
            {
                dictionary dict;

                dict.add("name", p.name());
                dict.add("type", "boundary");
                dict.add("T", eye);
                dict.add("boundaryOffset", offset);

                partPatches_.append(partPatch::New(*this, dict));

                goto found;
            }
        }

        FatalErrorInFunction
            << "Could not find the patch to which face " << facei
            << " of processor " << Pstream::myProcNo()
            << " belongs to."
            << endl << abort(FatalError);

        found:;
    }
}

void mesh::setCommTags()
{
    const label n = Pstream::nProcs();
    const label m = Pstream::myProcNo();

    List<List<dictionary>> patchDicts(n, List<dictionary>(0));

    forAll(partPatches_, patchi)
    {
        const partPatch& patch = partPatches_[patchi];

        if
        (
            patch.type() == "parallel"
        ||  patch.type() == "periodic"
        )
        {
            patchDicts[m].append
            (
                dictionary(patch.dict())
            );
        }
    }

    Pstream::gatherList(patchDicts);

    labelList myTags;

    if (Pstream::master())
    {
        List<labelList> tags(n,labelList(0));

        forAll(tags, i)
        {
            tags[i].setSize(patchDicts[i].size());
            tags[i] = 1;
        }

        // The pairCount array can get pretty big for a large number of
        // processors, but allocating an n by n array is by far the easiest way
        // to track pairs. On 10.000 processors, the array will consume 381 MB,
        // which should be doable.

        int pairCount[n][n];

        for (label i = 0; i < n; i++)
            for (label j = 0; j < n; j++)
                pairCount[i][j] = 0;

        forAll(patchDicts, i)
        forAll(patchDicts[i], patchi)
        {
            const dictionary& dict0 = patchDicts[i][patchi];

            const labelVector source
            (
              - labelTensor(dict0.lookup("T")).T()
              & labelVector(dict0.lookup("boundaryOffset"))
            );

            const label j = readLabel(dict0.lookup("neighborProcNum"));

            // Avoid searching pairs twice

            if (i <= j)
            {
                bool found = false;

                forAll(patchDicts[j], patchj)
                if (i != j || patchi != patchj)
                {
                    const dictionary& dict1 = patchDicts[j][patchj];

                    const label k =
                        readLabel(dict1.lookup("neighborProcNum"));

                    const labelVector target =
                        dict1.lookup("boundaryOffset");

                    if (k == i && source == target)
                    {
                        pairCount[i][j]++;

                        tags[i][patchi] = pairCount[i][j];
                        tags[j][patchj] = pairCount[i][j];

                        found = true;
                        break;
                    }
                }

                if (!found)
                {
                    FatalErrorInFunction
                        << "Could not find neighbor of patch "
                        << word(dict0.lookup("name"))
                        << " of processor " << j << endl
                        << abort(FatalError);
                }
            }
        }

        // Send

        forAll(patchDicts, i)
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

    forAll(partPatches_, patchi)
    {
        partPatch& patch = partPatches_[patchi];

        if
        (
            patch.type() == "parallel"
        ||  patch.type() == "periodic"
        )
        {
            patch.dict().add("tag", myTags[c++]);
        }
    }
}

void mesh::generatePartPatches()
{
    generateBrickInternalPartPatches();
    generateBrickExternalPartPatches();
    generateBoundaryPartPatches();

    setCommTags();

    lowerPatchMaster_ = unitXYZ;
    upperPatchMaster_ = unitXYZ;

    forAll(partPatches(), i)
    {
        const partPatch& patch = partPatches_[i];

        if (patch.slave() && patch.boundaryOffsetDegree() == 1)
        {
            const label facei = faceNumber(patch.boundaryOffset());

            if (facei % 2 == 0)
            {
                lowerPatchMaster_[facei/2] = 0;
            }
            else
            {
                upperPatchMaster_[facei/2] = 0;
            }
        }
    }
}

void mesh::generatePartLevels()
{
    label l = 0;
    bool add = true;
    const partLevel* parent = nullptr;

    while (add)
    {
        this->append
        (
            new partLevel(*this, parent)
        );

        parent = this->operator()(l);

        const labelVector P(parent->N());
        const labelVector D
        (
            P.x() < 4 ? P.x() : P.x()/2,
            P.y() < 4 ? P.y() : P.y()/2,
            P.z() < 4 ? P.z() : P.z()/2
        );

        label nProcsCoarsen = (P != D);

        reduce(nProcsCoarsen, sumOp<label>());

        if (nProcsCoarsen != Pstream::nProcs())
        {
            add = false;
        }

        l++;
    }
}

mesh::mesh(const IOdictionary& dict)
:
    geometry(dict),
    PtrList<partLevel>(0),
    decomp_(decomposition::New(*this)),
    N_(decomp_->myPartN()),
    lowerPatchMaster_(unitXYZ),
    upperPatchMaster_(unitXYZ)
{
    generatePartPatches();
    generatePartLevels();
}

mesh::~mesh()
{}

}

}
