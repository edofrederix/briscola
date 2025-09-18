#include "boundary.H"
#include "mesh.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(boundary, 0);
defineRunTimeSelectionTable(boundary, dictionary);

void boundary::extend()
{
    for (int i = 0; i < 6; i++)
        if (offset_[i/2] == 0)
            extension_[i] = 1;
}

boundary::boundary(const mesh& msh, const dictionary& dict)
:
    dict_(dict),
    name_(dict.lookup("name")),
    offset_(dict.lookup("offset")),
    offsetDegree_
    (
        cmptSum(cmptMag(offset_))
    ),
    extension_(Zero)
{}

boundary::boundary
(
    const boundary& b
)
:
    dict_(b.dict_),
    name_(b.name_),
    offset_(b.offset_),
    offsetDegree_(b.offsetDegree_),
    extension_(Zero)
{}

boundary::~boundary()
{}

autoPtr<boundary> boundary::New
(
    const mesh& msh,
    const dictionary& dict
)
{
    const word boundaryType(dict.lookup("type"));

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(boundaryType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown boundary type " << boundaryType << endl
            << "Valid boundary types are" << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<boundary>(cstrIter()(msh, dict));
}

}

}
