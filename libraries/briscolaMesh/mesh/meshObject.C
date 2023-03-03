#include "meshObject.H"

#include "edge.H"
#include "face.H"
#include "brick.H"

namespace Foam
{

namespace briscola
{

template<class ParentType>
meshObject<ParentType>::meshObject(const ParentType& p, const label num)
:
    parent_(p),
    num_(num)
{}

template<class ParentType>
meshObject<ParentType>::~meshObject()
{}

template class meshObject<edge>;
template class meshObject<face>;
template class meshObject<brick>;
template class meshObject<geometry>;

}

}
