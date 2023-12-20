#include "meshLevel.H"
#include "meshField.H"
#include "fvMesh.H"

#include "boundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::allocate()
{
    listType::setSize(MeshType::numberOfDirections);

    forAll(*this, d)
    {
        listType::set
        (
            d,
            new meshDirection<Type,MeshType>
            (
                *this,
                fvMsh_,
                l_,
                d
            )
        );
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::setLevelPointers()
{
    forAll(*this, d)
        listType::operator[](d).mshLevelPtr_ = this;
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::transfer
(
    meshLevel<Type,MeshType>& L
)
{
    listType::transfer(L);

    mshFieldPtr_ = L.mshFieldPtr_;
    L.mshFieldPtr_ = nullptr;

    setLevelPointers();
}

// Main constructor

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    meshField<Type,MeshType>& field,
    const fvMesh& fvMsh,
    const label l
)
:
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(fvMsh),
    l_(l),
    mshFieldPtr_(&field)
{
    allocate();
}

// Copy constructors

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const meshLevel<Type,MeshType>& L
)
:
    PtrList<meshDirection<Type,MeshType>>(L),
    refCount(),
    fvMsh_(L.fvMsh_),
    l_(L.levelNum()),
    mshFieldPtr_(nullptr)
{
    setLevelPointers();
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const meshLevel<Type,MeshType>& L,
    const zero&
)
:
    PtrList<meshDirection<Type,MeshType>>(L, Zero),
    refCount(),
    fvMsh_(L.fvMsh_),
    l_(L.levelNum()),
    mshFieldPtr_(nullptr)
{
    setLevelPointers();
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const meshLevel<Type,MeshType>& L,
    const Type& v
)
:
    PtrList<meshDirection<Type,MeshType>>(L, v),
    refCount(),
    fvMsh_(L.fvMsh_),
    l_(L.levelNum()),
    mshFieldPtr_(nullptr)
{
    setLevelPointers();
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const meshLevel<Type,MeshType>& L,
    const List<Type>& v
)
:
    PtrList<meshDirection<Type,MeshType>>(L, v),
    refCount(),
    fvMsh_(L.fvMsh_),
    l_(L.levelNum()),
    mshFieldPtr_(nullptr)
{
    setLevelPointers();
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const tmp<meshLevel<Type,MeshType>>& tL
)
:
    PtrList<meshDirection<Type,MeshType>>
    (
        const_cast<meshLevel<Type,MeshType>&>(tL()),
        tL.isTmp()
    ),
    refCount(),
    fvMsh_(tL->fvMsh_),
    l_(tL->levelNum()),
    mshFieldPtr_(nullptr)
{
    setLevelPointers();
    tL.clear();
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const tmp<meshLevel<Type,MeshType>>& tL,
    const zero&
)
:
    PtrList<meshDirection<Type,MeshType>>
    (
        const_cast<meshLevel<Type,MeshType>&>(tL()),
        tL.isTmp()
    ),
    refCount(),
    fvMsh_(tL->fvMsh_),
    l_(tL->levelNum()),
    mshFieldPtr_(nullptr)
{
    setLevelPointers();
    *this = Zero;
    tL.clear();
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const tmp<meshLevel<Type,MeshType>>& tL,
    const Type& v
)
:
    PtrList<meshDirection<Type,MeshType>>
    (
        const_cast<meshLevel<Type,MeshType>&>(tL()),
        tL.isTmp()
    ),
    refCount(),
    fvMsh_(tL->fvMsh_),
    l_(tL->levelNum()),
    mshFieldPtr_(nullptr)
{
    setLevelPointers();
    *this = v;
    tL.clear();
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const tmp<meshLevel<Type,MeshType>>& tL,
    const List<Type>& v
)
:
    PtrList<meshDirection<Type,MeshType>>
    (
        const_cast<meshLevel<Type,MeshType>&>(tL()),
        tL.isTmp()
    ),
    refCount(),
    fvMsh_(tL->fvMsh_),
    l_(tL->levelNum()),
    mshFieldPtr_(nullptr)
{
    setLevelPointers();
    *this = v;
    tL.clear();
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const fvMesh& fvMsh,
    const label l
)
:
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(fvMsh),
    l_(l),
    mshFieldPtr_(nullptr)
{
    allocate();
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const fvMesh& fvMsh,
    const label l,
    const zero&
)
:
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(fvMsh),
    l_(l),
    mshFieldPtr_(nullptr)
{
    allocate();
    *this = Zero;
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const fvMesh& fvMsh,
    const label l,
    const Type& v
)
:
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(fvMsh),
    l_(l),
    mshFieldPtr_(nullptr)
{
    allocate();
    *this = v;
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const fvMesh& fvMsh,
    const label l,
    const List<Type>& v
)
:
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(fvMsh),
    l_(l),
    mshFieldPtr_(nullptr)
{
    allocate();
    *this = v;
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::~meshLevel()
{}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctUnsetBoundaryConditions()
{
    forAll(*this, d)
    {
        meshDirection<Type,MeshType>& field = listType::operator[](d);

        forAll(fvMsh_.msh().emptyPatchOffsets(), i)
        {
            const labelVector bo(fvMsh_.msh().emptyPatchOffsets()[i]);

            const labelVector S(fvMsh_.template S<MeshType>(l_,d,bo));
            const labelVector E(fvMsh_.template E<MeshType>(l_,d,bo));

            labelVector ijk;

            for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
            for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
            for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
            {
                field(ijk+bo) = field(ijk);
            }
        }
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        // Correct unset boundary conditions first

        this->correctUnsetBoundaryConditions();

        // Next, correct all boundary conditions contained by this part

        const label nReq = Pstream::nRequests();

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            mshFieldPtr_->boundaryConditions()[i].prepare(l_);
        }

        if (Pstream::parRun())
        {
            Pstream::waitRequests(nReq);
        }

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            mshFieldPtr_->boundaryConditions()[i].evaluate(l_);
        }
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctParallelBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        // Correct all parallel boundary conditions

        const label nReq = Pstream::nRequests();

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            if (bc.baseType() == PARALLELBC)
            {
                bc.prepare(l_);
            }
        }

        if (Pstream::parRun())
        {
            Pstream::waitRequests(nReq);
        }

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            if (bc.baseType() == PARALLELBC)
            {
                bc.evaluate(l_);
            }
        }
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctPeriodicBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        // Correct all periodic boundary conditions

        const label nReq = Pstream::nRequests();

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            if (bc.baseType() == PERIODICBC)
            {
                bc.prepare(l_);
            }
        }

        if (Pstream::parRun())
        {
            Pstream::waitRequests(nReq);
        }

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            if (bc.baseType() == PERIODICBC)
            {
                bc.evaluate(l_);
            }
        }
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctCommBoundaryConditions()
{
    this->correctParallelBoundaryConditions();
    this->correctPeriodicBoundaryConditions();
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctNonCommBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        // Correct unset boundary conditions first

        this->correctUnsetBoundaryConditions();

        // Correct all non-communicating boundary conditions

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            if (bc.baseType() != PARALLELBC && bc.baseType() != PERIODICBC)
            {
                bc.prepare(l_);
            }
        }

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            if (bc.baseType() != PARALLELBC && bc.baseType() != PERIODICBC)
            {
                bc.evaluate(l_);
            }
        }
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator=(const meshLevel<Type,MeshType>& L)
{
    forAll(*this, d)
        listType::operator[](d) = L[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator=
(
    const tmp<meshLevel<Type,MeshType>>& tL
)
{
    if (tL.isTmp())
    {
        meshLevel<Type,MeshType>& L =
            const_cast<meshLevel<Type,MeshType>&>(tL());

        transfer(L);
    }
    else
    {
        *this = tL();
    }

    tL.clear();
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator=(const Type& v)
{
    forAll(*this, d)
        listType::operator[](d) = v;
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator=(const List<Type>& v)
{
    forAll(*this, d)
        listType::operator[](d) = v[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator=(const zero)
{
    forAll(*this, d)
        listType::operator[](d) = Zero;
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator+=(const meshLevel<Type,MeshType>& L)
{
    forAll(*this, d)
        listType::operator[](d) += L[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator+=
(
    const tmp<meshLevel<Type,MeshType>>& tL
)
{
    *this += tL();
    tL.clear();
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator-=(const meshLevel<Type,MeshType>& L)
{
    forAll(*this, d)
        listType::operator[](d) -= L[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator-=
(
    const tmp<meshLevel<Type,MeshType>>& tL
)
{
    *this -= tL();
    tL.clear();
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator*=(const meshLevel<scalar,MeshType>& L)
{
    forAll(*this, d)
        listType::operator[](d) *= L[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator*=
(
    const tmp<meshLevel<scalar,MeshType>>& tL
)
{
    *this *= tL();
    tL.clear();
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator/=(const meshLevel<scalar,MeshType>& L)
{
    forAll(*this, d)
        listType::operator[](d) /= L[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator/=
(
    const tmp<meshLevel<scalar,MeshType>>& tL
)
{
    *this /= tL();
    tL.clear();
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator+=(const Type& v)
{
    forAll(*this, d)
        listType::operator[](d) += v;
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator+=(const List<Type>& v)
{
    forAll(*this, d)
        listType::operator[](d) += v[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator-=(const Type& v)
{
    forAll(*this, d)
        listType::operator[](d) -= v;
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator-=(const List<Type>& v)
{
    forAll(*this, d)
        listType::operator[](d) -= v[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator*=(const scalar& v)
{
    forAll(*this, d)
        listType::operator[](d) *= v;
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator*=(const List<scalar>& v)
{
    forAll(*this, d)
        listType::operator[](d) *= v[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator/=(const scalar& v)
{
    forAll(*this, d)
        listType::operator[](d) /= v;
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator/=(const List<scalar>& v)
{
    forAll(*this, d)
        listType::operator[](d) /= v[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator=
(
    const meshLevel<Type2,MeshType>& L
)
{
    forAll(*this, d)
        listType::operator[](d) = L[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator=
(
    const tmp<meshLevel<Type2,MeshType>>& tL
)
{
    *this = tL();
    tL.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator+=
(
    const meshLevel<Type2,MeshType>& L
)
{
    forAll(*this, d)
        listType::operator[](d) += L[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator+=
(
    const tmp<meshLevel<Type2,MeshType>>& tL
)
{
    *this += tL();
    tL.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator-=
(
    const meshLevel<Type2,MeshType>& L
)
{
    forAll(*this, d)
        listType::operator[](d) -= L[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator-=
(
    const tmp<meshLevel<Type2,MeshType>>& tL
)
{
    *this -= tL();
    tL.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator*=
(
    const meshLevel<Type2,MeshType>& L
)
{
    forAll(*this, d)
        listType::operator[](d) *= L[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator*=
(
    const tmp<meshLevel<Type2,MeshType>>& tL
)
{
    *this *= tL();
    tL.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator/=
(
    const meshLevel<Type2,MeshType>& L
)
{
    forAll(*this, d)
        listType::operator[](d) /= L[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator/=
(
    const tmp<meshLevel<Type2,MeshType>>& tL
)
{
    *this /= tL();
    tL.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator=(const Type2& v)
{
    forAll(*this, d)
        listType::operator[](d) = v;
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator=(const List<Type2>& v)
{
    forAll(*this, d)
        listType::operator[](d) = v[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator+=(const Type2& v)
{
    forAll(*this, d)
        listType::operator[](d) += v;
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator+=(const List<Type2>& v)
{
    forAll(*this, d)
        listType::operator[](d) += v[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator-=(const Type2& v)
{
    forAll(*this, d)
        listType::operator[](d) -= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator-=(const List<Type2>& v)
{
    forAll(*this, d)
        listType::operator[](d) -= v[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator*=(const Type2& v)
{
    forAll(*this, d)
        listType::operator[](d) *= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator*=(const List<Type2>& v)
{
    forAll(*this, d)
        listType::operator[](d) *= v[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator/=(const Type2& v)
{
    forAll(*this, d)
        listType::operator[](d) /= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator/=(const List<Type2>& v)
{
    forAll(*this, d)
        listType::operator[](d) /= v[d];
}

}

}

}

#include "meshLevelFunctions.C"
#include "meshLevelStencilFunctions.C"
