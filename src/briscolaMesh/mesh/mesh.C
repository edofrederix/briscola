#include "mesh.H"

#include "unstructuredMesh.H"
#include "structuredMesh.H"
#include "rectilinearMesh.H"
#include "uniformMesh.H"
#include "periodicPatch.H"
#include "patchBoundary.H"
#include "emptyBoundary.H"
#include "parallelBoundary.H"
#include "periodicBoundary.H"

#include "PstreamGlobals.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(mesh, 0);
defineRunTimeSelectionTable(mesh, dictionary);

void mesh::addBoundary
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

    boundaries_.append(boundary::New(*this, dict));
}

void mesh::generateBrickInternalBoundaries()
{
    // Add parallel brick-internal boundaries

    const labelVector myBrickDecomp(decomp().myBrickDecomp());
    const labelVector myBrickPart(decomp().myBrickPart());

    if (cmptProduct(myBrickDecomp) > 1)
    {
        // Try all boundary offset vectors, including faces, edges and vertices.
        // Even if the bricks are unstructured, within a brick we can add edge
        // and vertex boundaries anyway.

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

            // Add to boundaries

            addBoundary
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

void mesh::generateBrickExternalBoundaries()
{
    const brickTopology& topo = topology();

    forAll(bricks(), bricki)
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

void mesh::generateLinkBoundaries(const brickLink& link)
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

        addBoundary
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
            vector::zero,
          - (T.T() & bo),
            slice.procNum1()
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
                    const brickDecompositionSlice& slice2 = slices[map(ijk2)];

                    addBoundary
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
                        vector::zero,
                      - (T.T() & bo2),
                        slice2.procNum1()
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
                    const brickDecompositionSlice& slice2 = slices[map(ijk2)];

                    addBoundary
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
                        vector::zero,
                      - (T.T() & bo2),
                        slice2.procNum1()
                    );
                }
            }
        }
    }
}

void mesh::generatePatchBoundaries()
{
    // Loop over all patches and find the corresponding brick face. If not
    // found, the faces of this part are not part of the patch and we set a fake
    // patch boundary associated with the patch, with a zero boundary offset.
    // This is needed to enable global parallel operations on patch boundaries.

    const label brickNum = decomp().myBrickNum();
    const brick& b = bricks()[brickNum];

    forAll(patches(), patchi)
    if (!patches()[patchi].castable<periodicPatch>())
    {
        const patch& p = patches()[patchi];

        bool found = false;

        // Multiple faces can belong to the same boundary, so don't break the
        // loop after finding a matching face.

        forAll(b.faces(), facei)
        {
            const face& fi = b.faces()[facei];

            // Skip faces that are already assigned to a boundary

            bool assigned = false;
            forAll(boundaries_, i)
                if (boundaries_[i].offset() == faceOffsets[facei])
                    assigned = true;

            if (assigned)
                continue;

            forAll(p.facePtrs(), facej)
            {
                const face& fj = *p.facePtrs()[facej];

                if (&fi == &fj)
                {
                    addBoundary
                    (
                        p.name(),
                        patches()[patchi].type(),
                        eye,
                        faceOffsets[facei],
                        p.dict().lookupOrDefault("offset",vector::zero)
                    );

                    found = true;
                }
            }
        }

        if (!found)
        {
            addBoundary
            (
                p.name(),
                patches()[patchi].type(),
                eye,
                zeroXYZ,
                p.dict().lookupOrDefault("offset",vector::zero)
            );
        }
    }

    // By now, all faces must have one and only one boundary associated with
    // them. Check.

    forAll(b.faces(), facei)
    {
        label count = 0;
        forAll(boundaries_, j)
            if (faceOffsets[facei] == boundaries_[j].offset())
                count++;

        if (count < 1)
        {
            FatalErrorInFunction
                << "Could not find the patch to which face " << facei
                << " belongs." << endl << abort(FatalError);
        }

        if (count > 1)
        {
            FatalErrorInFunction
                << "Face " << facei << " is associated with "
                << count << " boundaries." << endl << abort(FatalError);
        }
    }
}

void mesh::setCommTags()
{
    const label n = Pstream::nProcs();
    const label m = Pstream::myProcNo();

    List<List<dictionary>> bDicts(n, List<dictionary>(0));

    forAll(boundaries_, bi)
    {
        const boundary& b = boundaries_[bi];

        if (b.castable<parallelBoundary>())
        {
            bDicts[m].append
            (
                dictionary(b.dict())
            );
        }
    }

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

    forAll(boundaries_, bi)
    {
        boundary& b = boundaries_[bi];

        if (b.castable<parallelBoundary>())
            b.cast<parallelBoundary>().setTag(myTags[c++]);
    }
}

void mesh::generateBoundaries()
{
    generateBrickInternalBoundaries();
    generateBrickExternalBoundaries();
    generatePatchBoundaries();

    reorderBoundaries();

    setCommTags();
    setBoundaryMask();
}

void mesh::setBoundaryMask()
{
    boundaryMask_.setSize(3,3,3);
    pBoundaryMask_.setSize(3,3,3);

    boundaryMask_ = Zero;
    pBoundaryMask_ = Zero;

    // Set boundary mask

    forAll(boundaries_, i)
        if (boundaries_[i].offsetDegree() > 0)
            boundaryMask_(boundaries_[i].offset()+unitXYZ) = 1;

    // Set parallel/periodic boundary mask

    forAll(boundaries_, i)
        if (boundaries_[i].offsetDegree() > 0)
            if (boundaries_[i].castable<parallelBoundary>())
                pBoundaryMask_(boundaries_[i].offset()+unitXYZ) = 1;
}

void mesh::reorderBoundaries()
{
    labelList oldToNew(boundaries_.size());

    int j = 0;

    // Patch boundaries

    forAll(boundaries_, i)
        if (boundaries_[i].castable<patchBoundary>())
            oldToNew[i] = j++;

    // Parallel/periodic boundaries

    for (int d = 1; d <= 3; d++)
        forAll(boundaries_, i)
            if
            (
                boundaries_[i].castable<parallelBoundary>()
             && boundaries_[i].offsetDegree() == d
            )
                oldToNew[i] = j++;

    // Empty boundaries

    forAll(boundaries_, i)
        if (boundaries_[i].castable<emptyBoundary>())
            oldToNew[i] = j++;

    if (j != boundaries_.size())
        FatalErrorInFunction
            << "Could not order boundaries"
            << endl << abort(FatalError);

    boundaries_.reorder(oldToNew);
}

void mesh::generateLevels()
{
    label l = 0;
    bool add = true;
    const part* parent = nullptr;

    while (add)
    {
        this->append
        (
            new part(*this, parent)
        );

        parent = this->operator()(l);

        const labelVector P(parent->N());
        const labelVector C(this->coarsen(P));

        label nProcsCoarsen = (P != C);

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
    FastPtrList<part>(0),
    decomp_(decomposition::New(*this)),
    structured_(topology().structured()),
    rectilinear_(Zero),
    uniform_(Zero)
{
    // Add boundaries and levels

    generateBoundaries();
    generateLevels();

    // Set structured mesh properties

    if (structured_)
    {
        const part& p = this->operator[](0);

        for (int d = 0; d < 3; d++)
        {
            // Mesh is rectilinear in a direction if all parts are rectilinear
            // in that direction

            rectilinear_[d] =
                returnReduce(p.rectilinear()[d], minOp<label>());

            // Mesh is uniform in a direction if all parts are uniform in that
            // direction

            uniform_[d] =
                returnReduce(p.uniform()[d], minOp<label>());
        }
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
    FastPtrList<part>(msh),
    decomp_(msh.decomp_),
    boundaries_(msh.boundaries_),
    boundaryMask_(msh.boundaryMask_),
    pBoundaryMask_(msh.pBoundaryMask_),
    structured_(msh.structured_),
    rectilinear_(msh.rectilinear_),
    uniform_(msh.uniform_),
    boundingBox_(msh.boundingBox_)
{}

mesh::mesh(autoPtr<mesh>& mshPtr)
:
    geometry(mshPtr()),
    FastPtrList<part>(mshPtr(), true),
    decomp_(mshPtr->decomp_, true),
    boundaries_(mshPtr->boundaries_, true),
    boundaryMask_(mshPtr->boundaryMask_),
    pBoundaryMask_(mshPtr->pBoundaryMask_),
    structured_(mshPtr->structured_),
    rectilinear_(mshPtr->rectilinear_),
    uniform_(mshPtr->uniform_),
    boundingBox_(mshPtr->boundingBox_)
{
    mshPtr.clear();
}

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
    const part& p = this->operator[](l);

    if
    (
        point.x() < p.boundingBox().left()
     || point.x() > p.boundingBox().right()
     || point.y() < p.boundingBox().bottom()
     || point.y() > p.boundingBox().top()
     || point.z() < p.boundingBox().aft()
     || point.z() > p.boundingBox().fore()
    )
    {
        // Not in bounding box of part

        return -unitXYZ;
    }
    else if (l == this->size()-1)
    {
        // Final level, linear search

        forAllBlock(p, i, j, k)
            if (p.points().pointInCell(point, i, j, k))
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
            // search if the mesh is not rectilinear. In this case it may be on
            // an edge.

            if (cmptProduct(p.rectilinear()) == 0)
            {
                forAllBlock(p, i, j, k)
                    if (p.points().pointInCell(point, i, j, k))
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
                if (p.points().pointInCell(point, i, j, k))
                    return labelVector(i,j,k);

            // Otherwise, search in the vicinity too.

            S = cmptMax(S-unitXYZ, zeroXYZ);
            E = cmptMin(E+unitXYZ, p.points().shape()-unitXYZ);

            for (int i = S.x(); i < E.x(); i++)
            for (int j = S.y(); j < E.y(); j++)
            for (int k = S.z(); k < E.z(); k++)
                if (p.points().pointInCell(point, i, j, k))
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

labelVector mesh::coarsen(const labelVector P) const
{
    if (!topology().structured() || !Pstream::parRun())
    {
        // Keep at least two cells in each direction on unstructured meshes, to
        // avoid heavy distortion

        return labelVector
        (
            P.x() <= 3 ? P.x() : P.x()/2,
            P.y() <= 3 ? P.y() : P.y()/2,
            P.z() <= 3 ? P.z() : P.z()/2
        );
    }
    else
    {
        // Refine until one or three cells in each direction

        labelVector Q
        (
            (P.x() == 1 || P.x() == 3) ? P.x() : P.x()/2,
            (P.y() == 1 || P.y() == 3) ? P.y() : P.y()/2,
            (P.z() == 1 || P.z() == 3) ? P.z() : P.z()/2
        );

        // Avoid having one cell in a periodic direction

        for (int d = 0; d < 3; d++)
        {
            if
            (
                Q[d] < 2
             && b(faceOffsets[2*d  ]).castable<periodicBoundary>()
             && b(faceOffsets[2*d+1]).castable<periodicBoundary>()
            )
            {
                Q[d] = 2;
            }
        }

        return Q;
    }
}

}

}
