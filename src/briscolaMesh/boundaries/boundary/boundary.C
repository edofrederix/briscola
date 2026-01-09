#include "boundary.H"
#include "level.H"

namespace Foam
{

namespace briscola
{

defineTypeNameAndDebug(boundary, 0);
defineRunTimeSelectionTable(boundary, dictionary);

boundary::boundary(const level& lvl, const dictionary& dict)
:
    lvl_(lvl),
    dict_(dict),
    name_(dict.lookup("name")),
    offset_(dict.lookup("offset")),
    offsetDegree_(cmptSum(cmptMag(offset_))),
    master_(true)
{}

boundary::boundary(const boundary& b)
:
    lvl_(b.lvl_),
    dict_(b.dict_),
    name_(b.name_),
    offset_(b.offset_),
    offsetDegree_(b.offsetDegree_),
    master_(b.master_)
{}

boundary::boundary(const boundary& b, const level& lvl)
:
    lvl_(lvl),
    dict_(b.dict_),
    name_(b.name_),
    offset_(b.offset_),
    offsetDegree_(b.offsetDegree_),
    master_(b.master_)
{}

boundary::~boundary()
{}

autoPtr<boundary> boundary::New
(
    const level& lvl,
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

    return autoPtr<boundary>(cstrIter()(lvl, dict));
}

}

}
