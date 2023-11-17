#include "mesh.H"

#include "unstructuredMesh.H"
#include "structuredMesh.H"
#include "rectilinearMesh.H"
#include "uniformMesh.H"

#include "boundaryPartPatch.H"
#include "emptyPartPatch.H"
#include "parallelPartPatch.H"
#include "periodicPartPatch.H"

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
                    patches()[patchi].type() == patch::EMPTY
                  ? "empty" : "boundary",
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
            patch.typeNum() == parallelPartPatch::typeNumber
         || patch.typeNum() == periodicPartPatch::typeNumber
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

        List<labelList> pairCount(n, labelList(n));

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
            patch.typeNum() == parallelPartPatch::typeNumber
         || patch.typeNum() == periodicPartPatch::typeNumber
        )
        {
            patch.dict().add("tag", myTags[c++]);
        }
    }
}

void mesh::setPatchExtension()
{
    edgePatchExtension_ = Zero;
    vertexPatchExtension_ = Zero;

    forAll(partPatches_, i)
    {
        const partPatch& patch = partPatches_[i];

        //  Add edges to extended faces, and vertices to extended edges

        for (int j = 0; j < 6; j++)
        if (patch.extension()[j] == 1)
        {
            const labelVector offset =
                patch.boundaryOffset() + faceOffsets[j];

            if (edgeNumber(offset) != -1)
            {
                edgePatchExtension_[edgeNumber(offset)] = 1;
            }
            else
            {
                vertexPatchExtension_[vertexNumber(offset)] = 1;
            }
        }

        // Also add vertices for face patches extended in two directions

        if (patch.boundaryOffsetDegree() == 1)
        {
            for (int j = 0; j < 5; j++)
            for (int k = j+1; k < 6; k++)
            if (patch.extension()[j] == 1 && patch.extension()[k] == 1)
            {
                const labelVector offset =
                    patch.boundaryOffset() + faceOffsets[j] + faceOffsets[k];

                if (vertexNumber(offset) != -1)
                    vertexPatchExtension_[vertexNumber(offset)] = 1;
            }
        }
    }
}

void mesh::setEmptyPatchOffsets()
{
    emptyPatchOffsets_.clear();

    for (int i = 0; i < 12; i++)
        if (edgePatchType()[i] == -1 && edgePatchExtension()[i] == 0)
            emptyPatchOffsets_.append(edgeOffsets[i]);

    for (int i = 0; i < 8; i++)
        if (vertexPatchType()[i] == -1 && vertexPatchExtension()[i] == 0)
            emptyPatchOffsets_.append(vertexOffsets[i]);
}

void mesh::setPatchLabels()
{
    // Initialize edge and vertex master label to -1 as their corresponding part
    // patches may not exist

    facePatchMasterPerProc_.resize
    (
        Pstream::nProcs(),
        pTraits<faceLabel>::one
    );

    edgePatchMasterPerProc_.resize
    (
        Pstream::nProcs(),
      - pTraits<edgeLabel>::one
    );

    vertexPatchMasterPerProc_.resize
    (
        Pstream::nProcs(),
      - pTraits<vertexLabel>::one
    );

    // Initialize edge and vertex patch types to -1 as their corresponding part
    // patches may not exist

    facePatchTypePerProc_.resize
    (
        Pstream::nProcs(),
        pTraits<faceLabel>::zero
    );

    edgePatchTypePerProc_.resize
    (
        Pstream::nProcs(),
      - pTraits<edgeLabel>::one
    );

    vertexPatchTypePerProc_.resize
    (
        Pstream::nProcs(),
      - pTraits<vertexLabel>::one
    );

    forAll(partPatches_, i)
    {
        const partPatch& patch = partPatches_[i];

        if (patch.boundaryOffsetDegree() == 1)
        {
            const label facei = faceNumber(patch.boundaryOffset());

            facePatchMasterPerProc_[Pstream::myProcNo()][facei] =
                patch.master();

            facePatchTypePerProc_[Pstream::myProcNo()][facei] =
                patch.typeNum();
        }
        else if (patch.boundaryOffsetDegree() == 2)
        {
            const label edgei = edgeNumber(patch.boundaryOffset());

            edgePatchMasterPerProc_[Pstream::myProcNo()][edgei] =
                patch.master();

            edgePatchTypePerProc_[Pstream::myProcNo()][edgei] =
                patch.typeNum();
        }
        else
        {
            const label vertexi = vertexNumber(patch.boundaryOffset());

            vertexPatchMasterPerProc_[Pstream::myProcNo()][vertexi] =
                patch.master();

            vertexPatchTypePerProc_[Pstream::myProcNo()][vertexi] =
                patch.typeNum();
        }
    }

    Pstream::gatherList(facePatchMasterPerProc_);
    Pstream::gatherList(edgePatchMasterPerProc_);
    Pstream::gatherList(vertexPatchMasterPerProc_);

    Pstream::gatherList(facePatchTypePerProc_);
    Pstream::gatherList(edgePatchTypePerProc_);
    Pstream::gatherList(vertexPatchTypePerProc_);

    Pstream::scatterList(facePatchMasterPerProc_);
    Pstream::scatterList(edgePatchMasterPerProc_);
    Pstream::scatterList(vertexPatchMasterPerProc_);

    Pstream::scatterList(facePatchTypePerProc_);
    Pstream::scatterList(edgePatchTypePerProc_);
    Pstream::scatterList(vertexPatchTypePerProc_);
}

void mesh::generatePartPatches()
{
    generateBrickInternalPartPatches();
    generateBrickExternalPartPatches();
    generateBoundaryPartPatches();

    setCommTags();
    setPatchLabels();
    setPatchExtension();
    setEmptyPatchOffsets();
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
    decomp_(decomposition::New(*this))
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

    // Determine bounding box

    const faceScalar bb(this->operator[](0).boundingBox());

    boundingBox_ =
        faceScalar
        (
            returnReduce(bb.left(), minOp<scalar>()),
            returnReduce(bb.right(), maxOp<scalar>()),
            returnReduce(bb.bottom(), minOp<scalar>()),
            returnReduce(bb.top(), maxOp<scalar>()),
            returnReduce(bb.aft(), minOp<scalar>()),
            returnReduce(bb.fore(), maxOp<scalar>())
        );
}

mesh::mesh(const mesh& msh)
:
    geometry(msh),
    PtrList<partLevel>(msh),
    decomp_(msh.decomp_),
    partPatches_(msh.partPatches_),
    facePatchMasterPerProc_(msh.facePatchMasterPerProc_),
    edgePatchMasterPerProc_(msh.edgePatchMasterPerProc_),
    vertexPatchMasterPerProc_(msh.vertexPatchMasterPerProc_),
    facePatchTypePerProc_(msh.facePatchTypePerProc_),
    edgePatchTypePerProc_(msh.edgePatchTypePerProc_),
    vertexPatchTypePerProc_(msh.vertexPatchTypePerProc_),
    edgePatchExtension_(msh.edgePatchExtension_),
    vertexPatchExtension_(msh.vertexPatchExtension_),
    emptyPatchOffsets_(msh.emptyPatchOffsets_),
    structured_(msh.structured_),
    rectilinear_(msh.rectilinear_),
    uniform_(msh.uniform_),
    boundingBox_(msh.boundingBox_)
{}

mesh::mesh(autoPtr<mesh>& mshPtr)
:
    geometry(mshPtr(), true),
    PtrList<partLevel>(mshPtr(), true),
    decomp_(mshPtr->decomp_, true),
    partPatches_(mshPtr->partPatches_, true),
    facePatchMasterPerProc_(mshPtr->facePatchMasterPerProc_),
    edgePatchMasterPerProc_(mshPtr->edgePatchMasterPerProc_),
    vertexPatchMasterPerProc_(mshPtr->vertexPatchMasterPerProc_),
    facePatchTypePerProc_(mshPtr->facePatchTypePerProc_),
    edgePatchTypePerProc_(mshPtr->edgePatchTypePerProc_),
    vertexPatchTypePerProc_(mshPtr->vertexPatchTypePerProc_),
    edgePatchExtension_(mshPtr->edgePatchExtension_),
    vertexPatchExtension_(mshPtr->vertexPatchExtension_),
    emptyPatchOffsets_(mshPtr->emptyPatchOffsets_),
    structured_(mshPtr->structured_),
    rectilinear_(mshPtr->rectilinear_),
    uniform_(mshPtr->uniform_),
    boundingBox_(mshPtr->boundingBox_)
{
    mshPtr.clear();
}

mesh::mesh(mesh& msh, bool reuse)
:
    geometry(msh, reuse),
    PtrList<partLevel>(msh, reuse),
    decomp_(msh.decomp_, reuse),
    partPatches_(msh.partPatches_, reuse),
    facePatchMasterPerProc_(msh.facePatchMasterPerProc_),
    edgePatchMasterPerProc_(msh.edgePatchMasterPerProc_),
    vertexPatchMasterPerProc_(msh.vertexPatchMasterPerProc_),
    facePatchTypePerProc_(msh.facePatchTypePerProc_),
    edgePatchTypePerProc_(msh.edgePatchTypePerProc_),
    vertexPatchTypePerProc_(msh.vertexPatchTypePerProc_),
    edgePatchExtension_(msh.edgePatchExtension_),
    vertexPatchExtension_(msh.vertexPatchExtension_),
    emptyPatchOffsets_(msh.emptyPatchOffsets_),
    structured_(msh.structured_),
    rectilinear_(msh.rectilinear_),
    uniform_(msh.uniform_),
    boundingBox_(msh.boundingBox_)
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

labelVector mesh::findCell(const vector& point, const label l) const
{
    const partLevel& lvl = this->operator[](l);

    if
    (
        point.x() < lvl.boundingBox().left()
     || point.x() > lvl.boundingBox().right()
     || point.y() < lvl.boundingBox().bottom()
     || point.y() > lvl.boundingBox().top()
     || point.z() < lvl.boundingBox().aft()
     || point.z() > lvl.boundingBox().fore()
    )
    {
        // Not in bounding box of part

        return -unitXYZ;
    }
    else if (l == this->size()-1)
    {
        // Final level, linear search

        forAllBlock(lvl, i, j, k)
            if (lvl.points().pointInCell(point, i, j, k))
                return labelVector(i,j,k);

        return -unitXYZ;
    }
    else
    {
        const labelVector coarse = this->findCell(point, l+1);
        const labelVector R = this->operator[](l+1).R();

        if (coarse == -unitXYZ)
        {
            // In bounding box, but not found on the coarse level. Try linear
            // search if the mesh is not rectlinear. In this case it may be on
            // an edge.

            if (cmptProduct(lvl.rectilinear()) == 0)
            {
                forAllBlock(lvl, i, j, k)
                    if (lvl.points().pointInCell(point, i, j, k))
                        return labelVector(i,j,k);
            }

            return -unitXYZ;
        }
        else
        {
            // Found on coarse level. Search enclosed cells only.

            labelVector S = cmptMultiply(coarse,R);
            labelVector E = cmptMultiply(coarse,R) + R;

            for (int i = S.x(); i < E.x(); i++)
            for (int j = S.y(); j < E.y(); j++)
            for (int k = S.z(); k < E.z(); k++)
                if (lvl.points().pointInCell(point, i, j, k))
                    return labelVector(i,j,k);

            // Otherwise, search in the vicinity too.

            S = cmptMax(S-unitXYZ, zeroXYZ);
            E = cmptMin(E+unitXYZ, lvl.points().shape()-unitXYZ);

            for (int i = S.x(); i < E.x(); i++)
            for (int j = S.y(); j < E.y(); j++)
            for (int k = S.z(); k < E.z(); k++)
                if (lvl.points().pointInCell(point, i, j, k))
                    return labelVector(i,j,k);

            // Otherwise, there is something wrong

            FatalErrorInFunction
                << "Found point " << point << " on level " << l+1
                << " in cell " << coarse << " but could not locate"
                << " a cell that contains the same point on level " << l
                << endl << abort(FatalError);

            return -unitXYZ;
        }
    }
}

List<labelVector> mesh::findCells(const vectorList& points, const label l) const
{
    List<labelVector> res(points.size());

    forAll(points, i)
    {
        res[i] = this->findCell(points[i], l);
    }

    return res;
}

}

}
