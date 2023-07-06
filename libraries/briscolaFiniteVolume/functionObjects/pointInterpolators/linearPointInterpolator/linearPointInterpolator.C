#include "linearPointInterpolator.H"
#include "geometry.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(linearPointInterpolator, 0);

addToRunTimeSelectionTable
(
    pointInterpolator,
    linearPointInterpolator,
    dictionary
);

linearPointInterpolator::linearPointInterpolator
(
    const fvMesh& fvMsh,
    const vectorList& points
)
:
    pointInterpolator(fvMsh, points),
    weights_(points.size())
{
    const meshDirection<vector,colocated>& cc =
        fvMsh_.template metrics<colocated>().cellCenters()[0][0];

    forAll(points_, i)
    {
        if (indices_[i] != -unitXYZ)
        {
            const label ii = indices_[i].x();
            const label jj = indices_[i].y();
            const label kk = indices_[i].z();

            const vector u(cellCoordinates_[i]);

            const label di = u.x() < 0.5 ? -1 : 1;
            const label dj = u.y() < 0.5 ? -1 : 1;
            const label dk = u.z() < 0.5 ? -1 : 1;

            const vertexVector vertices
            (
                cc(ii,   jj,   kk   ),
                cc(ii+di,jj,   kk   ),
                cc(ii,   jj+dj,kk   ),
                cc(ii+di,jj+dj,kk   ),
                cc(ii,   jj,   kk+dk),
                cc(ii+di,jj,   kk+dk),
                cc(ii,   jj+dj,kk+dk),
                cc(ii+di,jj+dj,kk+dk)
            );

            const vector v(interpolationWeights(points_[i],vertices));

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

linearPointInterpolator::linearPointInterpolator
(
    const linearPointInterpolator& interp
)
:
    pointInterpolator(interp),
    weights_(interp.weights_)
{}

#define INTERPFUNC(TYPE)                                                    \
List<TYPE> linearPointInterpolator::operator()                              \
(                                                                           \
    const meshField<TYPE,colocated>& field                                  \
)                                                                           \
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

template<class Type>
List<Type> linearPointInterpolator::interp
(
    const meshField<Type,colocated>& field
)
{
    const meshDirection<Type,colocated>& f = field[0][0];

    List<Type> values(points_.size());

    forAll(points_, i)
    {
        if (indices_[i] != -unitXYZ)
        {
            const label ii = indices_[i].x();
            const label jj = indices_[i].y();
            const label kk = indices_[i].z();

            const vector u(cellCoordinates_[i]);

            const label di = u.x() < 0.5 ? -1 : 1;
            const label dj = u.y() < 0.5 ? -1 : 1;
            const label dk = u.z() < 0.5 ? -1 : 1;

            values[i] =
                weights_[i].lba()*f(ii,   jj,   kk   )
              + weights_[i].rba()*f(ii+di,jj,   kk   )
              + weights_[i].lta()*f(ii,   jj+dj,kk   )
              + weights_[i].rta()*f(ii+di,jj+dj,kk   )
              + weights_[i].lbf()*f(ii,   jj,   kk+dk)
              + weights_[i].rbf()*f(ii+di,jj,   kk+dk)
              + weights_[i].ltf()*f(ii,   jj+dj,kk+dk)
              + weights_[i].rtf()*f(ii+di,jj+dj,kk+dk);
        }
        else
        {
            values[i] = Zero;
        }
    }

    return values;
}

}

}

}
