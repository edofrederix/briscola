#include "linearPointInterpolator.H"
#include "geometry.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTemplateTypeNameAndDebugWithName
(
    linearPointInterpolator<colocated>,
    "linear",
    0
);

defineTemplateTypeNameAndDebugWithName
(
    linearPointInterpolator<staggered>,
    "linear",
    0
);

addTemplatedToRunTimeSelectionTable
(
    pointInterpolator,
    linearPointInterpolator,
    colocated,
    dictionary
);

addTemplatedToRunTimeSelectionTable
(
    pointInterpolator,
    linearPointInterpolator,
    staggered,
    dictionary
);

template<class MeshType>
linearPointInterpolator<MeshType>::linearPointInterpolator
(
    const fvMesh& fvMsh,
    const vectorList& points,
    const bool global,
    const label l,
    const label d
)
:
    pointInterpolator<MeshType>(fvMsh, points, global, l, d),
    weights_(this->points_.size())
{
    const meshField<vector,MeshType>& cc =
        this->fvMsh_.template metrics<MeshType>().cellCenters();

    forAll(this->points_, i)
    {
        if (this->indices_[i] != -unitXYZ)
        {
            const label ii = this->indices_[i].x();
            const label jj = this->indices_[i].y();
            const label kk = this->indices_[i].z();

            const vector u(this->cellCoordinates_[i]);

            const label di = u.x() < 0.5 ? -1 : 1;
            const label dj = u.y() < 0.5 ? -1 : 1;
            const label dk = u.z() < 0.5 ? -1 : 1;

            const vertexVector vertices
            (
                cc(l,d,ii,   jj,   kk   ),
                cc(l,d,ii+di,jj,   kk   ),
                cc(l,d,ii,   jj+dj,kk   ),
                cc(l,d,ii+di,jj+dj,kk   ),
                cc(l,d,ii,   jj,   kk+dk),
                cc(l,d,ii+di,jj,   kk+dk),
                cc(l,d,ii,   jj+dj,kk+dk),
                cc(l,d,ii+di,jj+dj,kk+dk)
            );

            const vector v(interpolationWeights(this->points_[i],vertices));

            weights_[i] =
                vertexScalar
                (
                    (1-v.x())*(1-v.y())*(1-v.z()),
                    (  v.x())*(1-v.y())*(1-v.z()),
                    (1-v.x())*(  v.y())*(1-v.z()),
                    (  v.x())*(  v.y())*(1-v.z()),
                    (1-v.x())*(1-v.y())*(  v.z()),
                    (  v.x())*(1-v.y())*(  v.z()),
                    (1-v.x())*(  v.y())*(  v.z()),
                    (  v.x())*(  v.y())*(  v.z())
                );
        }
        else
        {
            weights_[i] = -vertexScalar::one;
        }
    }
}

template<class MeshType>
linearPointInterpolator<MeshType>::linearPointInterpolator
(
    const linearPointInterpolator<MeshType>& interp
)
:
    pointInterpolator<MeshType>(interp),
    weights_(interp.weights_)
{}

template<class MeshType>
linearPointInterpolator<MeshType>::~linearPointInterpolator()
{}

#define INTERPFUNC(TYPE)                                                    \
template<class MeshType>                                                    \
List<TYPE> linearPointInterpolator<MeshType>::operator()                    \
(                                                                           \
    const meshField<TYPE,MeshType>& field                                   \
) const                                                                     \
{                                                                           \
    return this->interp(field);                                             \
}

INTERPFUNC(scalar);
INTERPFUNC(vector);
INTERPFUNC(tensor);
INTERPFUNC(sphericalTensor);
INTERPFUNC(symmTensor);
INTERPFUNC(diagTensor);
INTERPFUNC(faceScalar);
INTERPFUNC(faceVector);

#undef INTERPFUNC

template<class MeshType>
template<class Type>
List<Type> linearPointInterpolator<MeshType>::interp
(
    const meshField<Type,MeshType>& field
) const
{
    List<Type> values(this->points_.size());

    forAll(this->points_, i)
    {
        if (this->indices_[i] != -unitXYZ)
        {
            const label ii = this->indices_[i].x();
            const label jj = this->indices_[i].y();
            const label kk = this->indices_[i].z();

            const vector u(this->cellCoordinates_[i]);

            const label di = u.x() < 0.5 ? -1 : 1;
            const label dj = u.y() < 0.5 ? -1 : 1;
            const label dk = u.z() < 0.5 ? -1 : 1;

            values[i] =
                weights_[i].lba()*field(this->l_,this->d_,ii,   jj,   kk   )
              + weights_[i].rba()*field(this->l_,this->d_,ii+di,jj,   kk   )
              + weights_[i].lta()*field(this->l_,this->d_,ii,   jj+dj,kk   )
              + weights_[i].rta()*field(this->l_,this->d_,ii+di,jj+dj,kk   )
              + weights_[i].lbf()*field(this->l_,this->d_,ii,   jj,   kk+dk)
              + weights_[i].rbf()*field(this->l_,this->d_,ii+di,jj,   kk+dk)
              + weights_[i].ltf()*field(this->l_,this->d_,ii,   jj+dj,kk+dk)
              + weights_[i].rtf()*field(this->l_,this->d_,ii+di,jj+dj,kk+dk);
        }
        else
        {
            values[i] = Zero;
        }
    }

    return this->combine(values);
}

// Instantiate

template class linearPointInterpolator<colocated>;
template class linearPointInterpolator<staggered>;

}

}

}
