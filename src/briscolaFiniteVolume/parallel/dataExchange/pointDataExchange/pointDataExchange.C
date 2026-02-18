#include "pointDataExchange.H"
#include "patchBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTemplateTypeNameAndDebug(pointDataExchange<colocated>, 0);
defineTemplateTypeNameAndDebug(pointDataExchange<staggered>, 0);

template<class MeshType>
pointDataExchange<MeshType>::pointDataExchange
(
    const List<vector>& points,
    const fvMesh& fvMsh,
    const label l,
    const label d
)
:
    dataExchange<MeshType>(fvMsh,l,d)
{
    init(points);
}

template<class MeshType>
pointDataExchange<MeshType>::pointDataExchange
(
    const pointDataExchange& e
)
:
    dataExchange<MeshType>(e),
    interpPtr_(e.interpPtr_->clone())
{}

template<class MeshType>
pointDataExchange<MeshType>::~pointDataExchange()
{}

template<class MeshType>
void pointDataExchange<MeshType>::init(const List<vector>& points)
{
    interpPtr_.reset
    (
        pointInterpolator<MeshType>::New
        (
            this->fvMsh_,
            points,
            "linear",
            false,
            this->l_,
            this->d_
        ).ptr()
    );

    const vectorList missing(interpPtr_->missingPoints());

    if (missing.size())
        FatalErrorInFunction
            << "Not all interpolation points were found. Missing points are: "
            << endl << missing << endl
            << "at the " << MeshType::typeName << " mesh level "
            << interpPtr_->l() << " direction " << interpPtr_->d() << endl
            << abort(FatalError);
}

template<class MeshType>
template<class Type>
List<Type> pointDataExchange<MeshType>::dataFunc
(
    const meshField<Type,MeshType>& field
) const
{
    return interpPtr_->operator()(field);
}

// Instantiate class and data functions

template class pointDataExchange<colocated>;
template class pointDataExchange<staggered>;

#define INSTANTIATE(TYPE,MESHTYPE)                                             \
template List<TYPE> pointDataExchange<MESHTYPE>::dataFunc                      \
(                                                                              \
    const meshField<TYPE,MESHTYPE>&                                            \
) const;

INSTANTIATE(scalar,colocated)
INSTANTIATE(vector,colocated)
INSTANTIATE(tensor,colocated)
INSTANTIATE(sphericalTensor,colocated)
INSTANTIATE(symmTensor,colocated)
INSTANTIATE(diagTensor,colocated)
INSTANTIATE(faceScalar,colocated)
INSTANTIATE(faceVector,colocated)

INSTANTIATE(scalar,staggered)
INSTANTIATE(vector,staggered)
INSTANTIATE(tensor,staggered)
INSTANTIATE(sphericalTensor,staggered)
INSTANTIATE(symmTensor,staggered)
INSTANTIATE(diagTensor,staggered)
INSTANTIATE(faceScalar,staggered)
INSTANTIATE(faceVector,staggered)

#undef INSTANTIATE

}

}

}
