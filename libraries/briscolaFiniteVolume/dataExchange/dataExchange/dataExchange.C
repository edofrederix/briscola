#include "dataExchange.H"
#include "boundaryPartPatch.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTemplateTypeNameAndDebug(dataExchange<colocated>, 0);
defineTemplateTypeNameAndDebug(dataExchange<staggered>, 0);

template<class MeshType>
dataExchange<MeshType>::dataExchange
(
    const fvMesh& fvMsh,
    const label l,
    const label d
)
:
    fvMsh_(fvMsh),
    l_(l),
    d_(d)
{}

template<class MeshType>
dataExchange<MeshType>::dataExchange
(
    const dataExchange& e
)
:
    fvMsh_(e.fvMsh_),
    l_(e.l_),
    d_(e.d_)
{}

template<class MeshType>
dataExchange<MeshType>::~dataExchange()
{}

// Instantiate

template class dataExchange<colocated>;
template class dataExchange<staggered>;

}

}

}
