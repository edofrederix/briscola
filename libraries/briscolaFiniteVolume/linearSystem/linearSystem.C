#include "linearSystem.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::transferData
(
    linearSystem<SType,Type,MeshType>& sys
)
{
    A_.transferData(sys.A_);
    b_.transferData(sys.b_);

    singular_ = sys.singular_;
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    meshField<Type, MeshType>& x
)
:
    tmp<linearSystem<SType,Type,MeshType>>::refCount(),
    fvMsh_(x.fvMsh()),
    x_(x),
    A_
    (
        IOobject::groupName("A", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    singular_(false)
{}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    const linearSystem<SType,Type,MeshType>& sys
)
:
    tmp<linearSystem<SType,Type,MeshType>>::refCount(),
    fvMsh_(sys.fvMsh_),
    x_(sys.x_),
    A_
    (
        IOobject::groupName("A", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    singular_(sys.singular_)
{
    *this = sys;
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
:
    tmp<linearSystem<SType,Type,MeshType>>::refCount(),
    fvMsh_(tSys->fvMsh_),
    x_(tSys->x_),
    A_
    (
        IOobject::groupName("A", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    singular_(tSys->singular_)
{
    if (tSys.isTmp())
    {
        linearSystem<SType,Type,MeshType>& sys =
            const_cast<linearSystem<SType,Type,MeshType>&>(tSys());

        transferData(sys);
    }
    else
    {
        *this = tSys();
    }

    tSys.clear();
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    const linearSystem<SType,Type,MeshType>& sys,
    meshField<Type, MeshType>& x
)
:
    tmp<linearSystem<SType,Type,MeshType>>::refCount(),
    fvMsh_(sys.fvMsh_),
    x_(x),
    A_
    (
        IOobject::groupName("A", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    singular_(sys.singular_)
{
    *this = sys;
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    meshField<Type, MeshType>& x
)
:
    tmp<linearSystem<SType,Type,MeshType>>::refCount(),
    fvMsh_(tSys->fvMsh_),
    x_(x),
    A_
    (
        IOobject::groupName("A", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    singular_(tSys->singular_)
{
    if (tSys.isTmp())
    {
        linearSystem<SType,Type,MeshType>& sys =
            const_cast<linearSystem<SType,Type,MeshType>&>(tSys());

        transferData(sys);
    }
    else
    {
        *this = tSys();
    }

    tSys.clear();
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::~linearSystem()
{}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::residual
(
    meshField<Type, MeshType>& res
) const
{
    forAll(res, l)
        this->residual(res[l]);
}

template<class SType, class Type, class MeshType>
tmp<meshField<Type, MeshType>>
linearSystem<SType,Type,MeshType>::residual() const
{
    tmp<meshField<Type,MeshType>> tRes
    (
        new meshField<Type,MeshType>
        (
            "residual",
            fvMsh_
        )
    );

    this->residual(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::residual
(
    meshLevel<Type,MeshType>& res
) const
{
    forAll(res, d)
        this->residual(res[d]);
}

template<class SType, class Type, class MeshType>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::residual(const label l) const
{
    tmp<meshLevel<Type,MeshType>> tRes
    (
        new meshLevel<Type,MeshType>(fvMsh_,l)
    );

    this->residual(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::residual
(
    meshDirection<Type,MeshType>& res
) const
{
    const label l = res.levelNum();
    const label d = res.directionNum();

    Amul(res, this->A()[l][d], x_[l][d]);

    res *= -1.0;
    res += this->b()[l][d];
}

template<class SType, class Type, class MeshType>
tmp<meshDirection<Type,MeshType>>
linearSystem<SType,Type,MeshType>::residual
(
    const label l,
    const label d
) const
{
    tmp<meshDirection<Type,MeshType>> tRes
    (
        new meshDirection<Type,MeshType>(fvMsh_,l,d)
    );

    this->residual(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::correctBoundaryConditions()
{
    this->correctBoundarySources();
    this->correctBoundaryCoefficients();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::correctBoundaryConditions
(
    const label l
)
{
    this->correctBoundarySources(l);
    this->correctBoundaryCoefficients(l);
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::correctBoundarySources()
{
    forAll(x_, l)
        this->correctBoundarySources(l);
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::correctBoundarySources
(
    const label l
)
{
    const meshLevel<Type,MeshType>& x = x_[l];
    const meshLevel<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes()[l];

    meshLevel<Type,MeshType>& b = this->b()[l];

    forAll(x_.boundaryConditions(), boundaryi)
    {
        const boundaryCondition<Type,MeshType>& bc =
            x_.boundaryConditions()[boundaryi];

        if (!bc.eliminated())
        {
            const labelVector bo(bc.boundaryOffset());

            forAll(x, d)
            if (x[d].shifted(bo))
            {
                const meshDirection<Type,MeshType>& xd = x[d];
                const meshDirection<scalar,MeshType>& cvd = cv[d];

                meshDirection<Type,MeshType>& bd = b[d];

                const labelVector S(bd.boundaryStart(bo));
                const labelVector E(bd.boundaryEnd(bo));

                labelVector ijk;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    bd(ijk) = cvd(ijk)*xd(ijk);
                }
            }
        }
    }
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::correctBoundaryCoefficients()
{
    forAll(x_, l)
        this->correctBoundaryCoefficients(l);
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::correctBoundaryCoefficients
(
    const label l
)
{
    const meshLevel<Type,MeshType>& x = x_[l];
    const meshLevel<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes()[l];

    meshLevel<SType,MeshType>& A = this->A()[l];

    forAll(x_.boundaryConditions(), boundaryi)
    {
        const boundaryCondition<Type,MeshType>& bc =
            x_.boundaryConditions()[boundaryi];

        if (!bc.eliminated())
        {
            const labelVector bo(bc.boundaryOffset());

            forAll(x, d)
            if (x[d].shifted(bo))
            {
                const meshDirection<scalar,MeshType>& cvd = cv[d];

                meshDirection<SType,MeshType>& Ad = A[d];

                const labelVector S(Ad.boundaryStart(bo));
                const labelVector E(Ad.boundaryEnd(bo));

                labelVector ijk;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    Ad(ijk) = diagStencil(cvd(ijk));
                }
            }
        }
    }
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator=
(
    const linearSystem<SType,Type,MeshType>& sys
)
{
    this->A() = sys.A();
    this->b() = sys.b();
    singular_ = sys.singular();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator=
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    if (tSys.isTmp())
    {
        linearSystem<SType,Type,MeshType>& sys =
            const_cast<linearSystem<SType,Type,MeshType>&>(tSys());

        transferData(sys);
    }
    else
    {
        *this = tSys();
    }

    tSys.clear();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const linearSystem<SType,Type,MeshType>& sys
)
{
    this->A() += sys.A();
    this->b() += sys.b();
    singular_ = sys.singular() ? singular_ : false;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    *this += tSys();
    tSys.clear();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const linearSystem<SType,Type,MeshType>& sys
)
{
    this->A() -= sys.A();
    this->b() -= sys.b();
    singular_ = sys.singular() ? singular_ : false;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    *this -= tSys();
    tSys.clear();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator=
(
    const meshField<Type,MeshType>& field
)
{
    const meshField<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes();

    this->A() = diagStencil(0.0);
    this->b() = -cv*field;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator=
(
    const tmp<meshField<Type,MeshType>>& tField
)
{
    *this = tField();
    tField.clear();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const meshField<Type,MeshType>& field
)
{
    const meshField<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes();

    this->b() -= cv*field;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const tmp<meshField<Type,MeshType>>& tField
)
{
    *this += tField();
    tField.clear();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const meshField<Type,MeshType>& field
)
{
    const meshField<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes();

    this->b() += cv*field;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const tmp<meshField<Type,MeshType>>& tField
)
{
    *this -= tField();
    tField.clear();
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator=
(
    const linearSystem<SType2,Type,MeshType>& sys
)
{
    this->A() = sys.A();
    this->b() = sys.b();
    singular_ = sys.singular();
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator=
(
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys
)
{
    *this = tSys();
    tSys.clear();
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const linearSystem<SType2,Type,MeshType>& sys
)
{
    this->A() += sys.A();
    this->b() += sys.b();
    singular_ = sys.singular() ? singular_ : false;
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys
)
{
    *this += tSys();
    tSys.clear();
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const linearSystem<SType2,Type,MeshType>& sys
)
{
    this->A() -= sys.A();
    this->b() -= sys.b();
    singular_ = sys.singular() ? singular_ : false;
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys
)
{
    *this -= tSys();
    tSys.clear();
}

}

}

}
