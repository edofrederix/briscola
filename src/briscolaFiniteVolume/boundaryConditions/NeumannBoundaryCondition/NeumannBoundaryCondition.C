#include "NeumannBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshFields.H"
#include "faceFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

template<class Type>
void NeumannBoundaryCondition<Type,colocated>::eliminateGhosts
(
    linearSystem<stencil,Type,colocated>& sys,
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label fd = f/2;

    meshField<stencil,colocated>& A = sys.A();
    meshField<Type,colocated>& b = sys.b();

    const faceField<scalar,colocated>& delta = this->faceDeltas();

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        const scalar ghostCoeff = A(l,d,ijk)[f+1];

        A(l,d,ijk)[0] += ghostCoeff;
        A(l,d,ijk)[f+1] = Zero;

        if (l == 0)
        {
            b(l,d,ijk) -=
                ghostCoeff*this->boundaryGradients_[l](ijk-S)
              / delta[fd](l,d,upperFaceNeighbor(ijk,f));
        }
    }
}

template<class Type>
void NeumannBoundaryCondition<Type,colocated>::evaluate
(
    const label l,
    const label d
)
{
    const scalar H = l == 0;

    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label fd = f/2;

    meshDirection<Type,colocated>& field = this->mshField_[l][d];

    const faceField<scalar,colocated>& delta = this->faceDeltas();

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        field(ijk+bo) =
            field(ijk)
          + H*this->boundaryGradients_[l](ijk-S)
          / delta[fd](l,d,upperFaceNeighbor(ijk,f));
    }
}

// Staggered

template<class Type>
void NeumannBoundaryCondition<Type,staggered>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    const label item = l*3 + d;

    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label g(faceNumber(-bo));
    const label fd = f/2;

    meshField<stencil,staggered>& A = sys.A();
    meshField<Type,staggered>& b = sys.b();

    const faceField<scalar,staggered>& delta = this->faceDeltas();

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;

    if (fd == d)
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar ghostCoeff = A(l,d,ijk)[f+1];

            A(l,d,ijk)[g+1] += ghostCoeff;
            A(l,d,ijk)[f+1] = Zero;

            const scalar dx2 =
                1.0/delta[fd](l,d,upperFaceNeighbor(ijk,f))
              + 1.0/delta[fd](l,d,upperFaceNeighbor(ijk,g));

            if (l == 0)
                b(l,d,ijk) -=
                    ghostCoeff*dx2*this->boundaryGradients_[item](ijk-S);
        }
    }
    else
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar ghostCoeff = A(l,d,ijk)[f+1];

            A(l,d,ijk)[0] += ghostCoeff;
            A(l,d,ijk)[f+1] = Zero;

            if (l == 0)
            {
                b(l,d,ijk) -=
                    ghostCoeff*this->boundaryGradients_[item](ijk-S)
                  / delta[fd](l,d,upperFaceNeighbor(ijk,f));
            }
        }
    }
}

template<class Type>
void NeumannBoundaryCondition<Type,staggered>::evaluate
(
    const label l,
    const label d
)
{
    const label item = l*3 + d;

    const labelVector bo(this->offset());
    const label f(faceNumber(bo));
    const label g(faceNumber(-bo));
    const label fd = f/2;

    const scalar H = l == 0;

    meshDirection<Type,staggered>& field = this->mshField_[l][d];

    const staggeredScalarFaceField& delta = this->faceDeltas();

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;

    if (fd == d)
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            const scalar dx2 =
                1.0/delta[fd](l,d,upperFaceNeighbor(ijk,f))
              + 1.0/delta[fd](l,d,upperFaceNeighbor(ijk,g));

            field(ijk+bo) =
                field(ijk-bo) + H*dx2*this->boundaryGradients_[item](ijk-S);
        }
    }
    else
    {
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            field(ijk+bo) =
                field(ijk)
              + H*this->boundaryGradients_[item](ijk-S)
              / delta[fd](l,d,upperFaceNeighbor(ijk,f));
        }
    }
}

}

}

}

