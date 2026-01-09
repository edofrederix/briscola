#include "meshDirection.H"
#include "meshField.H"
#include "boundaryCondition.H"
#include "patchBoundary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::allocate()
{
    // Allocate a block of its level shape with mesh type padding and ghost cell
    // layers on both sides

    blockType::reAllocate
    (
        fvMsh_[l_].N()
      + MeshType::padding[d_]
      + 2*ghosts*unitXYZ
    );
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::transfer
(
    meshDirection<Type,MeshType>& D
)
{
    levelPtr_ = D.levelPtr_;

    l_ = D.l_;
    d_ = D.d_;

    blockType::transfer(D);

    D.levelPtr_ = nullptr;
}

// Main constructor

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    meshLevel<Type,MeshType>& level,
    const label d
)
:
    block<Type>(),
    fvMsh_(level.fvMsh()),
    l_(level.levelNum()),
    d_(d),
    I_(fvMsh_.I<MeshType>(l_,d)),
    levelPtr_(&level)
{
    allocate();
}

// Copy constructors

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const meshDirection<Type,MeshType>& D
)
:
    block<Type>(D),
    fvMsh_(D.fvMsh_),
    l_(D.l_),
    d_(D.d_),
    I_(D.I_),
    levelPtr_(nullptr)
{}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const meshDirection<Type,MeshType>& D,
    const zero&
)
:
    block<Type>(D,Zero),
    fvMsh_(D.fvMsh_),
    l_(D.l_),
    d_(D.d_),
    I_(D.I_),
    levelPtr_(nullptr)
{}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const meshDirection<Type,MeshType>& D,
    const Type& v
)
:
    block<Type>(D,v),
    fvMsh_(D.fvMsh_),
    l_(D.l_),
    d_(D.d_),
    I_(D.I_),
    levelPtr_(nullptr)
{}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const meshDirection<Type,MeshType>& D,
    const List<Type>& v
)
:
    block<Type>(D, v[D.d_]),
    fvMsh_(D.fvMsh_),
    l_(D.l_),
    d_(D.d_),
    I_(D.I_),
    levelPtr_(nullptr)
{}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const tmp<meshDirection<Type,MeshType>>& tD
)
:
    block<Type>
    (
        tD.isTmp(),
        const_cast<meshDirection<Type,MeshType>&>(tD())
    ),
    fvMsh_(tD->fvMsh_),
    l_(tD->l_),
    d_(tD->d_),
    I_(tD->I_),
    levelPtr_(nullptr)
{
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const tmp<meshDirection<Type,MeshType>>& tD,
    const zero&
)
:
    block<Type>
    (
        tD.isTmp(),
        const_cast<meshDirection<Type,MeshType>&>(tD()),
        Zero
    ),
    fvMsh_(tD->fvMsh_),
    l_(tD->l_),
    d_(tD->d_),
    I_(tD->I_),
    levelPtr_(nullptr)
{
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const tmp<meshDirection<Type,MeshType>>& tD,
    const Type& v
)
:
    block<Type>
    (
        tD.isTmp(),
        const_cast<meshDirection<Type,MeshType>&>(tD()),
        v
    ),
    fvMsh_(tD->fvMsh_),
    l_(tD->l_),
    d_(tD->d_),
    I_(tD->I_),
    levelPtr_(nullptr)
{
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const tmp<meshDirection<Type,MeshType>>& tD,
    const List<Type>& v
)
:
    block<Type>
    (
        tD.isTmp(),
        const_cast<meshDirection<Type,MeshType>&>(tD()),
        v[tD->d_]
    ),
    fvMsh_(tD->fvMsh_),
    l_(tD->l_),
    d_(tD->d_),
    I_(tD->I_),
    levelPtr_(nullptr)
{
    if (tD.isTmp())
        tD.clear();
}

// Standalone constructors

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const fvMesh& fvMsh,
    const label l,
    const label d
)
:
    block<Type>(),
    fvMsh_(fvMsh),
    l_(l),
    d_(d),
    I_(fvMsh.I<MeshType>(l,d)),
    levelPtr_(nullptr)
{
    allocate();
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const fvMesh& fvMsh,
    const label l,
    const label d,
    const zero&
)
:
    block<Type>(),
    fvMsh_(fvMsh),
    l_(l),
    d_(d),
    I_(fvMsh.I<MeshType>(l,d)),
    levelPtr_(nullptr)
{
    allocate();
    *this = Zero;
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const fvMesh& fvMsh,
    const label l,
    const label d,
    const Type& v
)
:
    block<Type>(),
    fvMsh_(fvMsh),
    l_(l),
    d_(d),
    I_(fvMsh.I<MeshType>(l,d)),
    levelPtr_(nullptr)
{
    allocate();
    *this = v;
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const fvMesh& fvMsh,
    const label l,
    const label d,
    const List<Type>& v
)
:
    block<Type>(),
    fvMsh_(fvMsh),
    l_(l),
    d_(d),
    I_(fvMsh.I<MeshType>(l,d)),
    levelPtr_(nullptr)
{
    allocate();
    *this = v[d];
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::~meshDirection()
{}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator=
(
    const meshDirection<Type,MeshType>& D
)
{
    this->B() = D.B();
}

template<class Type, class MeshType>
tmp<meshDirection<typename meshDirection<Type,MeshType>::cmptType,MeshType>>
meshDirection<Type,MeshType>::component
(
    const label dir
) const
{
    tmp<meshDirection<cmptType,MeshType>> tD =
        meshDirection<cmptType,MeshType>::New(this->fvMsh_, this->l_, this->d_);

    tD.ref().B() = this->B().component(dir);

    return tD;
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::replace
(
    const label dir,
    const cmptType& v
)
{
    this->B().replace(dir,v);
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::replace
(
    const label dir,
    const meshDirection<cmptType,MeshType>& D
)
{
    this->B().replace(dir, D.B());
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::replace
(
    const label dir,
    const tmp<meshDirection<cmptType,MeshType>>& tD
)
{
    this->replace(dir,tD());

    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::max(const Type& v)
{
    this->B().max(v);
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::min(const Type& v)
{
    this->B().min(v);
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator=
(
    const tmp<meshDirection<Type,MeshType>>& tD
)
{
    if (tD.isTmp())
    {
        meshDirection<Type,MeshType>& D =
            const_cast<meshDirection<Type,MeshType>&>(tD());

        transfer(D);
    }
    else
    {
        *this = tD();
    }

    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator=(const Type& v)
{
    this->B() = v;
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator=(const zero)
{
    this->B() = Zero;
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator+=
(
    const meshDirection<Type,MeshType>& D
)
{
    this->B() += D.B();
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator+=
(
    const tmp<meshDirection<Type,MeshType>>& tD
)
{
    *this += tD();
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator-=
(
    const meshDirection<Type,MeshType>& D
)
{
    this->B() -= D.B();
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator-=
(
    const tmp<meshDirection<Type,MeshType>>& tD
)
{
    *this -= tD();
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator*=
(
    const meshDirection<scalar,MeshType>& D
)
{
    this->B() *= D.B();
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator*=
(
    const tmp<meshDirection<scalar,MeshType>>& tD
)
{
    *this *= tD();
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator/=
(
    const meshDirection<scalar,MeshType>& D
)
{
    this->B() /= D.B();
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator/=
(
    const tmp<meshDirection<scalar,MeshType>>& tD
)
{
    *this /= tD();
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator+=(const Type& v)
{
    this->B() += v;
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator-=(const Type& v)
{
    this->B() -= v;
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator*=(const scalar& v)
{
    this->B() *= v;
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::operator/=(const scalar& v)
{
    this->B() /= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator=
(
    const meshDirection<Type2,MeshType>& D
)
{
    this->B() = D.B();
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator=
(
    const tmp<meshDirection<Type2,MeshType>>& tD
)
{
    *this = tD();
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator+=
(
    const meshDirection<Type2,MeshType>& D
)
{
    this->B() += D.B();
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator+=
(
    const tmp<meshDirection<Type2,MeshType>>& tD
)
{
    *this += tD();
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator-=
(
    const meshDirection<Type2,MeshType>& D
)
{
    this->B() -= D.B();
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator-=
(
    const tmp<meshDirection<Type2,MeshType>>& tD
)
{
    *this -= tD();
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator*=
(
    const meshDirection<Type2,MeshType>& D
)
{
    this->B() *= D.B();
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator*=
(
    const tmp<meshDirection<Type2,MeshType>>& tD
)
{
    *this *= tD();
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator/=
(
    const meshDirection<Type2,MeshType>& D
)
{
    this->B() /= D.B();
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator/=
(
    const tmp<meshDirection<Type2,MeshType>>& tD
)
{
    *this /= tD();
    if (tD.isTmp())
        tD.clear();
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator=(const Type2& v)
{
    this->B() = v;
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator+=(const Type2& v)
{
    this->B() += v;
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator-=(const Type2& v)
{
    this->B() -= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator*=(const Type2& v)
{
    this->B() *= v;
}

template<class Type, class MeshType>
template<class Type2>
void meshDirection<Type,MeshType>::operator/=(const Type2& v)
{
    this->B() /= v;
}

}

}

}

#include "meshDirectionFunctions.C"
#include "meshDirectionStencilFunctions.C"
