#include "geometry.H"
#include "periodicPatch.H"

namespace Foam
{

namespace briscola
{

void geometry::checkPatchConsistency() const
{
    // Check for duplicate faces in patches

    for (label patchi = 0; patchi < patches_.size()-1; patchi++)
    {
        const List<const face*>& fsi = patches_[patchi].facePtrs();

        for (label facei = 0; facei < fsi.size(); facei++)
        {
            for (label patchj = patchi+1; patchj < patches_.size(); patchj++)
            {
                const List<const face*>& fsj = patches_[patchj].facePtrs();

                for (label facej = 0; facej < fsj.size(); facej++)
                {
                    if (fsi[facei] == fsj[facej])
                    {
                        FatalErrorInFunction
                            << "Found duplicate " << *fsi[facei]
                            << " in patches " << patches_[patchi].name()
                            << " and " << patches_[patchj].name()
                            << exit(FatalError);
                    }
                }
            }
        }
    }
}

void geometry::createBricks()
{
    bricks_.clear();
    bricks_.setSize(brickData().size());

    forAll(brickData(), i)
    {
        bricks_.set(i, new brick(*this, i, brickData()[i]));
    }
}

void geometry::createBrickTopology()
{
    topology_.clear();
    topology_.set(new brickTopology(*this));
}

void geometry::createPatches()
{
    patches_.clear();

    forAllConstIter(patchTable, patchData(), iter)
        addPatch(patch(*this, patches_.size(), iter.key(), iter()));
}

void geometry::createDefaultPatch()
{
    // Add all faces which are not part of a patch nor connected to another
    // brick to the default patch

    List<labelList> defaultFaces;

    for (label bricki = 0; bricki < bricks_.size(); bricki++)
    {
        const brick& b = bricks_[bricki];

        forAll(b.faces(), facei)
        {
            const face& f = b.faces()[facei];

            bool found = false;

            if (topology_->links()[bricki].faceLinks().set(facei))
            {
                found = true;
            }
            else
            {
                forAll(patches_, patchi)
                {
                    if (patches_[patchi].hasFace(f))
                    {
                        found = true;
                    }
                }
            }

            if (!found)
            {
                defaultFaces.append(f.vertexList());
            }
        }
    }

    if (defaultFaces.size() > 0)
    {
        dictionary dict;
        dict.add("type", patch::typeName);
        dict.add("faces", defaultFaces);

        addPatch(patch(*this, patches_.size(), "default", dict));
    }
}

void geometry::addPatch(const patch& p)
{
    // A patch may not contain two faces of the same brick. So, in some cases,
    // the original patch needs to be split up in multiple patches. Once a brick
    // is detected that contains two or more faces, faces are separated from
    // each other based on their face numbers in the brick.

    // Collect face and brick numbers

    labelList faceNums(p.facePtrs().size());
    labelList brickNums(p.facePtrs().size());

    forAll(p.facePtrs(), i)
    {
        faceNums[i] = p.facePtrs()[i]->num();
        brickNums[i] = p.facePtrs()[i]->parentBrick().num();
    }

    // Check number of unique brick numbers

    labelList order;
    uniqueOrder(brickNums, order);

    // If the number of unique brick numbers is not equal to the number of
    // bricks then there are faces that share the same brick. In that case we
    // split up the patch into patches with faces having the same face number in
    // their respective brick (i.e., they have the same orientation).

    List<List<labelList>> faces(6, List<labelList>(0));

    if (order.size() != brickNums.size())
    {
        forAll(p.facePtrs(), i)
            faces[faceNums[i]].append(p.facePtrs()[i]->vertexList());
    }
    else
    {
        // Just copy to the first entry
        forAll(p.facePtrs(), i)
            faces[0].append(p.facePtrs()[i]->vertexList());
    }

    // Count number of sub-patches

    label nSubPatches = 0;

    forAll(faces, i)
        if (faces[i].size())
            nSubPatches++;

    // Add (sub-)patches

    label j = 0;
    forAll(faces, i)
    {
        if (faces[i].size())
        {
            dictionary dict(p.dict());
            dict.set("faces", faces[i]);

            const word name
            (
                nSubPatches > 1
              ? IOobject::groupName(p.name(), Foam::name(j++))
              : p.name()
            );

            patches_.append
            (
                patch::New
                (
                    *this,
                    patches_.size(),
                    name,
                    dict
                )
            );
        }
    }
}

void geometry::createPatchPairs()
{
    patchPairs_.clear();

    label i = 0;

    for (label p0 = 0; p0 < patches_.size(); p0++)
    {
        if (patches_[p0].castable<periodicPatch>())
        {
            const periodicPatch& patch0 =
                patches_[p0].cast<periodicPatch>();

            const word neighbor0(patch0.neighbor());

            bool found = false;

            for (label p1 = 0; p1 < patches_.size(); p1++)
            {
                if (patches_[p1].name() == neighbor0)
                {
                    found = true;

                    const periodicPatch& patch1 =
                        patches_[p1].cast<periodicPatch>();

                    const word neighbor1(patch1.neighbor());

                    if (neighbor1 != patch0.name())
                    {
                        FatalErrorInFunction
                            << "Patch " << patch0.name()
                            << "'s periodic neighbor is "
                            << neighbor0 << ", but patch " << neighbor0
                            << "'s periodic "
                            << "neighbor is " << neighbor1 << endl
                            << exit(FatalError);
                    }

                    if (patch0.facePtrs().size() != patch1.facePtrs().size())
                    {
                        FatalErrorInFunction
                            << "Patch " << patch0.name()
                            << " forms a periodic patch pair "
                            << "with patch " << patch1.name()
                            << " but they have a different "
                            << "number of faces." << endl
                            << exit(FatalError);
                    }

                    bool foundPair = false;

                    forAll(patchPairs_, j)
                    {
                        const patchPair& pair = patchPairs_[j];

                        if
                        (
                            (pair.patch0() == patch0 && pair.patch1() == patch1)
                         || (pair.patch1() == patch0 && pair.patch0() == patch1)
                        )
                        {
                            foundPair = true;
                            break;
                        }
                    }

                    if (!foundPair)
                    {
                        patchPairs_.append
                        (
                            new patchPair(*this, i, patch0, patch1)
                        );
                        i++;
                    }

                    break;
                }
            }

            if (!found)
            {
                FatalErrorInFunction
                    << "Could not find neighboring patch named " << neighbor0
                    << " of periodic patch " << patch0.name() << endl
                    << exit(FatalError);
            }
        }
    }
}

void geometry::alignBricks()
{
    if (!topology_.valid())
    {
        FatalErrorInFunction
            << "Brick topology not set."
            << exit(FatalError);
    }

    if (!topology_->structured() || topology_->aligned() || bricks_.size() == 1)
    {
        return;
    }

    // Collect brick transformations first, because transforming them on the fly
    // invalidates the brick topology. Align bricks w.r.t. the orientation of
    // the largest brick. This means that all bricks will have a local
    // coordinate system that is aligned with that of the largest brick. This is
    // the orientation that is most likely to have most cells in the third local
    // direction, which gives the highest contiguous memory density.

    labelList sizes(bricks_.size());

    forAll(bricks_, bricki)
        sizes[bricki] = cmptProduct(bricks_[bricki].N());

    const label largest = findMax(sizes);

    List<labelTensor> transforms(bricks_.size());

    forAll(bricks_, bricki)
    {
        if (bricki != largest)
        {
            // Find path between the largest brick and this brick

            const labelList P = topology_->shortestFacePath(largest, bricki);

            labelTensor T = eye;

            for (label i = P.size()-2; i >= 0; i--)
            {
                const brickLinks& links = topology_->links()[P[i]];
                const brickFaceLink& link = links.getFaceLink(P[i+1]);

                T = link.T() & T;
            }

            transforms[bricki] = T;
        }
    }

    // Perform brick transform now that all transforms are known

    forAll(transforms, bricki)
        if (bricki != largest)
            if (transforms[bricki] != eye)
                bricks_[bricki].transform(transforms[bricki]);

    createPatches();
    createPatchPairs();
    checkPatchConsistency();
    createBrickTopology();

    // Check new topology

    if (!topology_->aligned())
        FatalErrorInFunction
            << "Brick alignment failed." << endl << abort(FatalError);
}

geometry::geometry(const IOdictionary& dict)
:
    meshData(dict),
    bricks_(),
    topology_(),
    patches_(),
    patchPairs_(),
    dataAlignment_
    (
        dict.lookupOrDefault<Switch>("dataAlignment", true)
    )
{
    createBricks();
    createPatches();
    createPatchPairs();
    checkPatchConsistency();
    createBrickTopology();
    alignBricks();
    createDefaultPatch();
    checkPatchConsistency();
}

geometry::geometry(const geometry& geo)
:
    meshData(geo),
    bricks_(geo.bricks_, *this),
    dataAlignment_(geo.dataAlignment_)
{
    createPatches();
    createPatchPairs();
    checkPatchConsistency();
    createBrickTopology();
    alignBricks();
    createDefaultPatch();
    checkPatchConsistency();
}

geometry::~geometry()
{}

}

}
