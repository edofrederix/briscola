#include "mesh.H"

#include "unstructuredMesh.H"
#include "structuredMesh.H"
#include "rectilinearMesh.H"
#include "uniformMesh.H"

#include "PstreamGlobals.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(mesh, 0);
defineRunTimeSelectionTable(mesh, dictionary);

void mesh::generateLevels()
{
    label l = 0;
    bool add = true;
    const level* parent = nullptr;

    while (add)
    {
        this->append(new level(*this, parent));

        parent = this->operator()(l);

        bool coarsenable = parent->coarsenable();

        reduce(coarsenable, andOp<bool>());

        if (!coarsenable)
        {
            // Not all levels are coarsenable. Check if we can coarsen the
            // decomposition. Empty levels are always coarsenable.

            coarsenable =
                parent->empty() || parent->decomp().coarsenable();

            reduce(coarsenable, andOp<bool>());

            // Stop if also the decomposition cannot be coarsened

            if (!coarsenable)
                add = false;
        }

        l++;
    }
}

mesh::mesh(const IOdictionary& dict)
:
    geometry(dict),
    FastPtrList<level>(0),
    structured_(topology().structured()),
    rectilinear_(Zero),
    uniform_(Zero)
{
    // Generate levels

    generateLevels();

    // Set structured mesh properties

    if (structured_)
    {
        const level& lvl = this->operator[](0);

        for (int d = 0; d < 3; d++)
        {
            // Mesh is rectilinear in a direction if all levels are rectilinear
            // in that direction

            rectilinear_[d] =
                returnReduce(lvl.rectilinear()[d], minOp<label>());

            // Mesh is uniform in a direction if all levels are uniform in that
            // direction

            uniform_[d] =
                returnReduce(lvl.uniform()[d], minOp<label>());
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
    FastPtrList<level>(msh, *this),
    structured_(msh.structured_),
    rectilinear_(msh.rectilinear_),
    uniform_(msh.uniform_),
    boundingBox_(msh.boundingBox_)
{}

mesh::mesh(autoPtr<mesh>& mshPtr)
:
    geometry(mshPtr()),
    FastPtrList<level>(mshPtr(), *this),
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
    const level& lvl = this->operator[](l);

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
        // Not in bounding box of level

        return -unitXYZ;
    }
    else if (lvl.leaf() || lvl.decomp().aggParent())
    {
        // Linear search on leaf or agglomerate parent levels

        forAllBlock(lvl, i, j, k)
            if (lvl.points().pointInCell(point, i, j, k))
                return labelVector(i,j,k);

        return -unitXYZ;
    }
    else
    {
        const labelVector coarse = this->findCell(point, l+1);
        const labelVector R = lvl.child().R();

        if (coarse == -unitXYZ)
        {
            // In bounding box, but not found on the coarse level. Try linear
            // search if the mesh is not rectilinear. In this case it may be on
            // an edge.

            if (cmptProduct(lvl.rectilinear()) == 0)
                forAllBlock(lvl, i, j, k)
                    if (lvl.points().pointInCell(point, i, j, k))
                        return labelVector(i,j,k);

            return -unitXYZ;
        }
        else
        {
            // Found on coarse level. Search enclosed cells only.

            labelVector S = cmptMax(cmptMultiply(coarse,R), zeroXYZ);
            labelVector E = cmptMin(S + R, lvl.N());

            for (int i = S.x(); i < E.x(); i++)
            for (int j = S.y(); j < E.y(); j++)
            for (int k = S.z(); k < E.z(); k++)
                if (lvl.points().pointInCell(point, i, j, k))
                    return labelVector(i,j,k);

            // Otherwise, search in the vicinity too.

            S = cmptMax(S - unitXYZ, zeroXYZ);
            E = cmptMin(E + unitXYZ, lvl.N());

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
