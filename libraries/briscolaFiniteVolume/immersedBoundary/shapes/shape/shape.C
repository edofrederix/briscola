#include "shape.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Constructor

shape::shape
(
    const dictionary& dict,
    bool inverted
)
:
    inverted_(inverted)
{}

autoPtr<shape> shape::New
(
    const dictionary& dict,
    bool inverted
)
{
    const word shapeType(dict.lookup("type"));

    typename dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(shapeType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown shape "
            << shapeType << nl << nl
            << "Valid solvers are:" << nl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<shape>(cstrIter()(dict, inverted));
}

// Destructor

shape::~shape()
{}

} // end namespace fv

} // end namespace briscola

} // end namespace Foam
