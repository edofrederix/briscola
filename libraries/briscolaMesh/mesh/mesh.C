#include "mesh.H"

#include "unstructuredMesh.H"
#include "structuredMesh.H"
#include "rectilinearMesh.H"
#include "uniformMesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(mesh, 0);
defineRunTimeSelectionTable(mesh, dictionary);

void mesh::addPartPatch
(
    const word name,
    const word type,
    const labelTensor T,
    const labelVector boundaryOffset,
    const labelVector neighborOffset,
    const label neighborProcNum
)
{
    dictionary dict;

    dict.add("name", name);
    dict.add("type", type);
    dict.add("T", T);
    dict.add("boundaryOffset", boundaryOffset);

    if (neighborOffset != zeroXYZ)
        dict.add("neighborOffset", neighborOffset);

    if (neighborProcNum >= 0)
        dict.add("neighborProcNum", neighborProcNum);

    partPatches_.append
    (
        partPatch::New(*this, dict)
    );
}
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

            addPartPatch
            (
                word
                (
                    "parallel-" + Foam::name(Pstream::myProcNo())
                  + "to" + Foam::name(neighborProcNum)
                ),
                "parallel",
                eye,
                bo,
              - bo,
                neighborProcNum
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

    const labelBlock& map = interface.map();
    const PtrList<brickDecompositionSlice>& slices = interface.slices();

    const labelVector bo = link.offset();
    const labelTensor T = link.T();

    // Only add slices if they belong to this processor

    forAllBlock(map, i, j, k)
    if (slices[map(i,j,k)].procNum0() == Pstream::myProcNo())
    {
        const labelVector ijk(i,j,k);
        const brickDecompositionSlice& slice = slices[map(ijk)];

        addPartPatch
        (
            word
            (
                word(link.periodic() ? "periodic" : "parallel")
              + "-" + Foam::name(slice.procNum0())
              + "to" + Foam::name(slice.procNum1())
            ),
            link.periodic() ? "periodic" : "parallel",
            T,
            bo,
          - (T.T() & bo),
            slice.procNum1()
        );

        // Add internal edges and vertices of face part patches

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
                    const brickDecompositionSlice& slice2 = slices[map(ijk2)];

                    addPartPatch
                    (
                        word
                        (
                            word(link.periodic() ? "periodic" : "parallel")
                          + "-" + Foam::name(slice.procNum0())
                          + "to" + Foam::name(slice2.procNum1())
                        ),
                        link.periodic() ? "periodic" : "parallel",
                        T,
                        bo2,
                      - (T.T() & bo2),
                        slice2.procNum1()
                    );
                }
            }
        }

        // Add internal vertices of edge part patches

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
                    const brickDecompositionSlice& slice2 = slices[map(ijk2)];

                    addPartPatch
                    (
                        word
                        (
                            word(link.periodic() ? "periodic" : "parallel")
                          + "-" + Foam::name(slice.procNum0())
                          + "to" + Foam::name(slice2.procNum1())
                        ),
                        link.periodic() ? "periodic" : "parallel",
                        T,
                        bo2,
                      - (T.T() & bo2),
                        slice2.procNum1()
                    );
                }
            }
        }
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
                addPartPatch
                (
                    p.name(),
                    "boundary",
                    eye,
                    offset
                );

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

    lowerMasterPatch_ = unitXYZ;
    upperMasterPatch_ = unitXYZ;

    forAll(partPatches(), i)
    {
        const partPatch& patch = partPatches_[i];

        if (patch.slave() && patch.boundaryOffsetDegree() == 1)
        {
            const label facei = faceNumber(patch.boundaryOffset());

            if (facei % 2 == 0)
            {
                lowerMasterPatch_[facei/2] = 0;
            }
            else
            {
                upperMasterPatch_[facei/2] = 0;
            }
        }
    }

    lowerBoundaryPatch_ = zeroXYZ;
    upperBoundaryPatch_ = zeroXYZ;

    forAll(partPatches(), i)
    {
        const partPatch& patch = partPatches_[i];

        if (patch.type() == "boundary" && patch.boundaryOffsetDegree() == 1)
        {
            const label facei = faceNumber(patch.boundaryOffset());

            if (facei % 2 == 0)
            {
                lowerBoundaryPatch_[facei/2] = 1;
            }
            else
            {
                upperBoundaryPatch_[facei/2] = 1;
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
    lowerMasterPatch_(unitXYZ),
    upperMasterPatch_(unitXYZ),
    lowerBoundaryPatch_(zeroXYZ),
    upperBoundaryPatch_(zeroXYZ)
{
    generatePartPatches();
    generatePartLevels();

    // Mesh is structured if the brick topology is too

    structured_ = topology().structured();

    if (structured_)
    {
        const partLevel& part = this->operator[](0);

        for (int d = 0; d < 3; d++)
        {
            // Mesh is rectilinear in a direction if all parts are rectilinear
            // in that direction

            rectilinear_[d] =
                returnReduce(part.rectilinear()[d], minOp<label>());

            // Mesh is uniform in a direction of all parts are uniform in that
            // direction

            uniform_[d] =
                returnReduce(part.uniform()[d], minOp<label>());
        }
    }
    else
    {
        rectilinear_ = zeroXYZ;
        uniform_ = zeroXYZ;
    }
}

mesh::mesh(const mesh& msh)
:
    geometry(msh),
    PtrList<partLevel>(msh),
    decomp_(msh.decomp_),
    partPatches_(msh.partPatches_),
    lowerMasterPatch_(msh.lowerMasterPatch_),
    upperMasterPatch_(msh.upperMasterPatch_),
    lowerBoundaryPatch_(msh.lowerBoundaryPatch_),
    upperBoundaryPatch_(msh.upperBoundaryPatch_),
    structured_(msh.structured_),
    rectilinear_(msh.rectilinear_),
    uniform_(msh.uniform_)
{}

mesh::mesh(autoPtr<mesh>& mshPtr)
:
    geometry(mshPtr(), true),
    PtrList<partLevel>(mshPtr(), true),
    decomp_(mshPtr->decomp_, true),
    partPatches_(mshPtr->partPatches_, true),
    lowerMasterPatch_(mshPtr->lowerMasterPatch_),
    upperMasterPatch_(mshPtr->upperMasterPatch_),
    lowerBoundaryPatch_(mshPtr->lowerBoundaryPatch_),
    upperBoundaryPatch_(mshPtr->upperBoundaryPatch_),
    structured_(mshPtr->structured_),
    rectilinear_(mshPtr->rectilinear_),
    uniform_(mshPtr->uniform_)
{
    mshPtr.clear();
}

mesh::mesh(mesh& msh, bool reuse)
:
    geometry(msh, reuse),
    PtrList<partLevel>(msh, reuse),
    decomp_(msh.decomp_, reuse),
    partPatches_(msh.partPatches_, reuse),
    lowerMasterPatch_(msh.lowerMasterPatch_),
    upperMasterPatch_(msh.upperMasterPatch_),
    lowerBoundaryPatch_(msh.lowerBoundaryPatch_),
    upperBoundaryPatch_(msh.upperBoundaryPatch_),
    structured_(msh.structured_),
    rectilinear_(msh.rectilinear_),
    uniform_(msh.uniform_)
{}

autoPtr<mesh> mesh::New(const IOdictionary& dict)
{
    // Construct the mesh as primitive type

    autoPtr<mesh> mshPtr(new mesh(dict));

    // Based on its properties, select the appropriate derived type

    word meshType;

    if (mshPtr->uniform() == unitXYZ)
    {
        meshType = uniformMesh::typeName;
    }
    else if (mshPtr->rectilinear() == unitXYZ)
    {
        meshType = rectilinearMesh::typeName;
    }
    else if (mshPtr->structured())
    {
        meshType = structuredMesh::typeName;
    }
    else
    {
        meshType = unstructuredMesh::typeName;
    }

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(meshType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown mesh type " << meshType << endl
            << ". Valid mesh types are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<mesh>(cstrIter()(mshPtr));
}

mesh::~mesh()
{}

}

}
