#include "slipBoundaryCondition.H"

#include "colocated.H"
#include "staggered.H"
#include "meshLevel.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Colocated

template<class Type>
void slipBoundaryCondition<Type,colocated>::evaluate
(
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    meshDirection<Type,colocated>& fd = this->mshField_[l][d];

    const meshDirection<faceVector,colocated>& fn = this->faceNormals()[l][d];

    const labelVector S(this->S(l,d));
    const labelVector E(this->E(l,d));

    labelVector ijk;
    for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
    for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
    for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
    {
        fd(ijk+bo) =
            fd(ijk) - 2.0*(fn(ijk)[faceNum] & fd(ijk)) * fn(ijk)[faceNum];
    }
}

// Staggered

template<class Type>
void slipBoundaryCondition<Type,staggered>::eliminateGhosts
(
    linearSystem<stencil,Type,staggered>& sys,
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    if (faceNum/2 == d)
    {
        meshField<stencil,staggered>& A = sys.A();
        meshField<Type,staggered>& b = sys.b();

        const meshField<scalar,staggered>& cv = this->cellVolumes();

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            A(l,d,ijk) = diagStencil(cv(l,d,ijk));

            if (l == 0)
                b(l,d,ijk) = Zero;
        }
    }
    else
    {
        NeumannBoundaryCondition<Type,staggered>::eliminateGhosts(sys,l,d);
    }
}

template<class Type>
void slipBoundaryCondition<Type,staggered>::evaluate
(
    const label l,
    const label d
)
{
    const labelVector bo(this->offset());
    const label faceNum(faceNumber(bo));

    if (faceNum/2 == d)
    {
        meshDirection<Type,staggered>& fd = this->mshField_[l][d];

        const labelVector S(this->S(l,d));
        const labelVector E(this->E(l,d));

        labelVector ijk;
        for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
        for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
        for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
        {
            fd(ijk+bo) = Zero;
        }
    }
    else
    {
        NeumannBoundaryCondition<Type,staggered>::evaluate(l,d);
    }
}

}

}

}

