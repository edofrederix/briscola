#include "meshLevel.H"
#include "meshField.H"
#include "fvMesh.H"

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
void meshLevel<Type,MeshType>::transferData
(
    meshLevel<Type,MeshType>& L
)
{
    listType::transfer(L);

    mshFieldPtr_ = L.mshFieldPtr_;
    L.mshFieldPtr_ = nullptr;

    forAll(*this, d)
        listType::operator[](d).mshLevelPtr_ = this;
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::updateActiveCells()
{
    forAll(*this, d)
    {
        listType::operator[](d).updateActiveCells();
    }
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
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(L.fvMsh_),
    l_(L.levelNum()),
    mshFieldPtr_(nullptr)
{
    allocate();
    *this = L;
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const meshLevel<Type,MeshType>& L,
    const zero&
)
:
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(L.fvMsh_),
    l_(L.levelNum()),
    mshFieldPtr_(nullptr)
{
    allocate();
    *this = Zero;
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const meshLevel<Type,MeshType>& L,
    const Type& v
)
:
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(L.fvMsh_),
    l_(L.levelNum()),
    mshFieldPtr_(nullptr)
{
    allocate();
    *this = v;
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const meshLevel<Type,MeshType>& L,
    const List<Type>& v
)
:
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(L.fvMsh_),
    l_(L.levelNum()),
    mshFieldPtr_(nullptr)
{
    allocate();
    *this = v;
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const tmp<meshLevel<Type,MeshType>>& tL
)
:
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(tL->fvMsh_),
    l_(tL->levelNum()),
    mshFieldPtr_(nullptr)
{
    if (tL.isTmp())
    {
        meshLevel<Type,MeshType>& L =
            const_cast<meshLevel<Type,MeshType>&>(tL());

        transferData(L);
    }
    else
    {
        allocate();
        *this = tL();
    }

    tL.clear();
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const tmp<meshLevel<Type,MeshType>>& tL,
    const zero&
)
:
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(tL->fvMsh_),
    l_(tL->levelNum()),
    mshFieldPtr_(nullptr)
{
    if (tL.isTmp())
    {
        meshLevel<Type,MeshType>& L =
            const_cast<meshLevel<Type,MeshType>&>(tL());

        transferData(L);
    }
    else
    {
        allocate();
    }

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
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(tL->fvMsh_),
    l_(tL->levelNum()),
    mshFieldPtr_(nullptr)
{
    if (tL.isTmp())
    {
        meshLevel<Type,MeshType>& L =
            const_cast<meshLevel<Type,MeshType>&>(tL());

        transferData(L);
    }
    else
    {
        allocate();
    }

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
    PtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(tL->fvMsh_),
    l_(tL->levelNum()),
    mshFieldPtr_(nullptr)
{
    if (tL.isTmp())
    {
        meshLevel<Type,MeshType>& L =
            const_cast<meshLevel<Type,MeshType>&>(tL());

        transferData(L);
    }
    else
    {
        allocate();
    }

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
void meshLevel<Type,MeshType>::initGhosts()
{
    initGhosts(pTraits<Type>::zero);
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::initGhosts(const Type& v)
{
    forAll(*this, d)
        listType::operator[](d).initGhosts(v);
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctBoundaryConditions
(
    const bool homogeneousBCs
)
{
    if (mshFieldPtr_)
    {
        // A call to correctBoundaryConditions() implies that boundary
        // conditions are needed for this level. Add them if not already done.

        if (mshFieldPtr_->boundaryConditions().size() == 0)
        {
            mshFieldPtr_->addBoundaryConditions();
        }

        // First, update all vertex and edge ghost cells as homogeneous Neumann.
        // This is needed because they might not be set by boundary conditions

        forAll(*this, d)
        {
            meshDirection<Type,MeshType>& field = listType::operator[](d);

            labelVector bo;

            for (bo.x() = -1; bo.x() <= 1; bo.x()++)
            for (bo.y() = -1; bo.y() <= 1; bo.y()++)
            for (bo.z() = -1; bo.z() <= 1; bo.z()++)
            if (cmptSum(cmptMag(bo)) > 1)
            {
                const labelVector S(field.boundaryStart(bo));
                const labelVector E(field.boundaryEnd(bo));

                labelVector ijk;

                for (ijk.x() = S.x(); ijk.x() < E.x(); ijk.x()++)
                for (ijk.y() = S.y(); ijk.y() < E.y(); ijk.y()++)
                for (ijk.z() = S.z(); ijk.z() < E.z(); ijk.z()++)
                {
                    field(ijk+bo) = field(ijk);
                }
            }
        }

        // Next, correct all boundary conditions contained by this part. This
        // may overwrite the edge and vertex Neumann BCs just set.

        const label nReq = Pstream::nRequests();

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            mshFieldPtr_->boundaryConditions()[i].initEvaluate(l_);
        }

        if (Pstream::parRun())
        {
            Pstream::waitRequests(nReq);
        }

        forAll(mshFieldPtr_->boundaryConditions(), i)
        {
            mshFieldPtr_->boundaryConditions()[i].evaluate(l_, homogeneousBCs);
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
void meshLevel<Type,MeshType>::operator=(const tmp<meshLevel<Type,MeshType>>& tL)
{
    if (tL.isTmp())
    {
        meshLevel<Type,MeshType>& L =
            const_cast<meshLevel<Type,MeshType>&>(tL());

        transferData(L);
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
void meshLevel<Type,MeshType>::operator+=(const tmp<meshLevel<Type,MeshType>>& tL)
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
void meshLevel<Type,MeshType>::operator-=(const tmp<meshLevel<Type,MeshType>>& tL)
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
void meshLevel<Type,MeshType>::operator*=(const tmp<meshLevel<scalar,MeshType>>& tL)
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
void meshLevel<Type,MeshType>::operator/=(const tmp<meshLevel<scalar,MeshType>>& tL)
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
void meshLevel<Type,MeshType>::operator=(const meshLevel<Type2,MeshType>& L)
{
    forAll(*this, d)
        listType::operator[](d) = L[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator=(const tmp<meshLevel<Type2,MeshType>>& tL)
{
    *this = tL();
    tL.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator+=(const meshLevel<Type2,MeshType>& L)
{
    forAll(*this, d)
        listType::operator[](d) += L[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator+=(const tmp<meshLevel<Type2,MeshType>>& tL)
{
    *this += tL();
    tL.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator-=(const meshLevel<Type2,MeshType>& L)
{
    forAll(*this, d)
        listType::operator[](d) -= L[d];
}

template<class Type, class MeshType>
template<class Type2>
void meshLevel<Type,MeshType>::operator-=(const tmp<meshLevel<Type2,MeshType>>& tL)
{
    *this -= tL();
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

}

}

}

#include "meshLevelFunctions.C"
#include "meshLevelStencilFunctions.C"
