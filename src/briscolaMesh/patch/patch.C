#include "patch.H"
#include "geometry.H"

namespace Foam
{

template<>
const char* NamedEnum<briscola::patch::patchType,4>::names[] =
{
    "patch",
    "periodic",
    "mapped",
    "empty"
};

const NamedEnum<briscola::patch::patchType,4> briscola::patch::patchTypeNames;

namespace briscola
{

void patch::checkConsistency() const
{
    for (label i = 0; i < facePtrs_.size()-1; i++)
    {
        for (label j = i+1; j < facePtrs_.size(); j++)
        {
            if (facePtrs_[i] == facePtrs_[j])
            {
                FatalErrorInFunction
                    << "Patch " << name_ << " has duplicate face"
                    << exit(FatalError);
            }
        }
    }
}

patch::patch
(
    const geometry& g,
    const label num,
    const word name,
    const dictionary& dict
)
:
    meshObject<geometry>(g, num),
    dictPtr_(&dict),
    name_(name),
    type_
    (
        patchTypeNames.read(dict.lookup("type"))
    ),
    facePtrs_()
{
    const List<labelList> vertexNumList(dict.lookup("faces"));

    if (vertexNumList.size() == 0)
    {
        FatalErrorInFunction
            << "No faces specified in patch " << name_
            << exit(FatalError);
    }

    facePtrs_.clear();
    facePtrs_.setSize(vertexNumList.size(), nullptr);

    forAll(vertexNumList, facei)
    {
        if (vertexNumList[facei].size() != 4)
        {
            FatalErrorInFunction
                << "Invalid face specified: " << vertexNumList[facei]
                << " for patch " << name
                << exit(FatalError);
        }

        const labelBlock vertexNums(2, 2, 1, vertexNumList[facei]);

        const face dummyFace
        (
            g.bricks()[0],
            0,
            vertexNums,
            unitXYZ
        );

        forAll(g.bricks(), bricki)
        {
            const brick& b = g.bricks()[bricki];

            forAll(b.faces(), facej)
            {
                if (dummyFace == b.faces()[facej])
                {
                    if (facePtrs_[facei] != nullptr)
                    {
                        FatalErrorInFunction
                            << "Cannot create patch from shared face"
                            << exit(FatalError);
                    }
                    else
                    {
                        facePtrs_[facei] = &b.faces()[facej];
                    }
                }
            }
        }

        if (facePtrs_[facei] == nullptr)
        {
            FatalErrorInFunction
                << "Could not find face " << vertexNumList[facei]
                << " of patch " << name
                << exit(FatalError);
        }
    }

    checkConsistency();
}

patch::patch
(
    const geometry& g,
    const label num,
    const word name,
    const List<const face*>& facePtrs
)
:
    meshObject<geometry>(g, num),
    dictPtr_(nullptr),
    name_(name),
    type_(PATCH),
    facePtrs_(facePtrs)
{
    checkConsistency();
}

patch::patch
(
    const patch& p
)
:
    meshObject<geometry>(p.parentGeometry(), p.num()),
    dictPtr_(&p.dict()),
    name_(p.name()),
    type_(p.type()),
    facePtrs_(p.facePtrs())
{}

patch::~patch()
{}

bool patch::operator==(const patch& p) const
{
    if (facePtrs_.size() != p.facePtrs().size())
    {
        return false;
    }

    forAll(facePtrs_, i)
    {
        const face& face0 = *facePtrs_[i];
        const face& face1 = *p.facePtrs()[i];

        if (face0 != face1)
        {
            return false;
        }
    }

    return true;
}

bool patch::operator!=(const patch& p) const
{
    return !(*this == p);
}

}

}
