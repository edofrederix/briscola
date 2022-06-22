#include "geometry.H"

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

    label i = 0;

    forAllConstIter(patchTable, patchData(), iter)
    {
        const dictionary& patchDict = iter();

        patches_.append(new patch(*this, i, iter.key(), patchDict));

        i++;
    }
}

void geometry::createDefaultPatch()
{
    // Add all faces which are not part of a patch nor connected to another
    // brick to the default patch

    List<const face*> defaultFaces;

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
                defaultFaces.append(&f);
            }
        }
    }

    if (defaultFaces.size() > 0)
    {
        patches_.append(new patch(*this, patches_.size(), "default", defaultFaces));
    }
}

void geometry::createPatchPairs()
{
    patchPairs_.clear();

    label i = 0;

    for (label p0 = 0; p0 < patches_.size(); p0++)
    {
        const patch& patch0 = patches_[p0];

        if (patch0.type() == patch::PERIODIC)
        {
            const word neighbor0(patch0.dict().lookup("neighbor"));

            bool found = false;

            for (label p1 = 0; p1 < patches_.size(); p1++)
            if (p1 != p0)
            {
                const patch& patch1 = patches_[p1];

                if (patch1.name() == neighbor0)
                {
                    found = true;

                    if (patch1.type() != patch::PERIODIC)
                    {
                        FatalErrorInFunction
                            << "Found periodic neighbor patch named " << patch1.name()
                            << " but it is not a periodic patch." << endl
                            << exit(FatalError);
                    }

                    const word neighbor1(patch1.dict().lookup("neighbor"));

                    if (neighbor1 != patch0.name())
                    {
                        FatalErrorInFunction
                            << "Patch " << patch0.name() << "'s periodic neighbor is "
                            << neighbor0 << ", but patch " << neighbor0 << "'s periodic "
                            << "neighbor is " << neighbor1 << endl
                            << exit(FatalError);
                    }

                    if (patch0.facePtrs().size() != patch1.facePtrs().size())
                    {
                        FatalErrorInFunction
                            << "Patch " << patch0.name() << " forms a periodic patch pair "
                            << "with patch " << patch1.name() << " but they have a different "
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
                        patchPairs_.append(new patchPair(*this, i, patch0, patch1));
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

geometry::geometry(const IOdictionary& dict)
:
    meshData(dict),
    dict_(dict),
    bricks_(),
    patches_(),
    patchPairs_()
{
    createBricks();
    createPatches();
    createPatchPairs();
    checkPatchConsistency();
    createBrickTopology();
    createDefaultPatch();
    checkPatchConsistency();
}

geometry::~geometry()
{}

}

}
