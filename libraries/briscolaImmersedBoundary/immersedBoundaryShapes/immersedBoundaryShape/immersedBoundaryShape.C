#include "immersedBoundaryShape.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

immersedBoundaryShape::immersedBoundaryShape
(
    const dictionary& dict,
    bool inverted
)
:
    inverted_(inverted)
{}

autoPtr<immersedBoundaryShape> immersedBoundaryShape::New
(
    const dictionary& dict,
    bool inverted
)
{
    const word immersedBoundaryShapeType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(immersedBoundaryShapeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown immersed boundary shape "
            << immersedBoundaryShapeType << nl << nl
            << "Valid shape types are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<immersedBoundaryShape>(cstrIter()(dict, inverted));
}

// Destructor

immersedBoundaryShape::~immersedBoundaryShape()
{}

}

}

}
