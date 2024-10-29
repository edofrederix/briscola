#include "meshLevel.H"
#include "meshField.H"
#include "fvMesh.H"

#include "boundaryCondition.H"
#include "immersedBoundaryCondition.H"

#include "domainBoundary.H"
#include "parallelBoundary.H"
#include "periodicBoundary.H"
#include "emptyBoundary.H"

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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
void meshLevel<Type,MeshType>::correctBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        // Immersed boundary conditions must be corrected first, so that cell
        // values near the boundaries can be copied to neighbors. By default,
        // only immersed boundary conditions on the first level are corrected.
        // If the immersed boundary conditions must be corrected on another
        // level, then the correctImmersedBoundaryConditions() must be called
        // explicitly on that level.

        if (l_ == 0)
            this->correctImmersedBoundaryConditions();

        this->correctUnsetBoundaryConditions();
        this->correctDomainBoundaryConditions();
        this->correctCommsBoundaryConditions();
        this->correctEmptyBoundaryConditions();
    }
}

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
void meshLevel<Type,MeshType>::correctDomainBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            const boundary& b = bc.mshBoundary();

            if (b.castable<domainBoundary>())
            {
                bc.evaluate(l_);
            }
        }
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctEmptyBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            const boundary& b = bc.mshBoundary();

            if (b.castable<emptyBoundary>())
            {
                bc.evaluate(l_);
            }
        }
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctCommsBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        #ifdef BOUNDARYEXCHANGE

        this->mshFieldPtr_->bExchangePtr_->correct(l_);

        #else

        const label nReq = Pstream::nRequests();

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            const boundary& b = bc.mshBoundary();

            if (b.castable<parallelBoundary>())
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

            const boundary& b = bc.mshBoundary();

            if (b.castable<parallelBoundary>())
            {
                bc.evaluate(l_);
            }
        }

        #endif
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctParallelBoundaryConditions()
{
    if (mshFieldPtr_ && Pstream::parRun())
    {
        this->addBoundaryConditions();

        #ifdef BOUNDARYEXCHANGE

        this->mshFieldPtr_->bExchangePtr_->correctParallel(l_);

        #else

        const label nReq = Pstream::nRequests();

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            const boundary& b = bc.mshBoundary();

            if
            (
                b.castable<parallelBoundary>()
            && !b.castable<periodicBoundary>()
            )
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

            const boundary& b = bc.mshBoundary();

            if
            (
                b.castable<parallelBoundary>()
            && !b.castable<periodicBoundary>()
            )
            {
                bc.evaluate(l_);
            }
        }

        #endif
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctPeriodicBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        #ifdef BOUNDARYEXCHANGE

        this->mshFieldPtr_->bExchangePtr_->correctPeriodic(l_);

        #else

        const label nReq = Pstream::nRequests();

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            const boundary& b = bc.mshBoundary();

            if (b.castable<periodicBoundary>())
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

            const boundary& b = bc.mshBoundary();

            if (b.castable<periodicBoundary>())
            {
                bc.evaluate(l_);
            }
        }

        #endif
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctNonEliminatedBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        // First non-eliminated domain boundaries

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            const boundary& b = bc.mshBoundary();

            if (!bc.eliminated() && b.castable<domainBoundary>())
            {
                bc.evaluate(l_);
            }
        }

        // Next the parallel/periodic boundaries which are non-eliminated

        this->correctCommsBoundaryConditions();

        // Finally the other non-eliminated boundaries

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            const boundary& b = bc.mshBoundary();

            if
            (
                !bc.eliminated()
             && !b.castable<domainBoundary>()
             && !b.castable<parallelBoundary>()
            )
            {
                bc.evaluate(l_);
            }
        }
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctEliminatedBoundaryConditions()
{
    if (mshFieldPtr_)
    {
        this->addBoundaryConditions();

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            boundaryCondition<Type,MeshType>& bc =
                mshFieldPtr_->boundaryConditions()[i];

            if (bc.eliminated())
            {
                bc.evaluate(l_);
            }
        }
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctImmersedBoundaryConditions()
{
    if (mshFieldPtr_ && fvMsh_.ibs<MeshType>().size())
    {
        this->addImmersedBoundaryConditions();

        forAll(mshFieldPtr_->immersedBoundaryConditions(), i)
        {
            immersedBoundaryCondition<Type,MeshType>& ibc =
                mshFieldPtr_->immersedBoundaryConditions()[i];

            ibc.evaluate(l_);
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

    if (tL.isTmp())
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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
    if (tL.isTmp())
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
