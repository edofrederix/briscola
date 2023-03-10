#include "meshDirection.H"
#include "meshField.H"
#include "boundaryCondition.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::setActiveCells()
{
    setActiveCells
    (
        fvMsh_.msh().lowerSlavePatch(),
        fvMsh_.msh().upperSlavePatch()
    );
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::setActiveCells
(
    const labelVector lowerSlave,
    const labelVector upperSlave
)
{
    const labelVector& padding = MeshType::padding[d_];

    S_ =
        briscola::cmptMultiply(padding, lowerSlave);

    E_ =
        this->B().N()
      - 2*unitXYZ
      - briscola::cmptMultiply(padding, upperSlave);

    N_ = E_ - S_;
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::updateActiveCells()
{
    const PtrList<boundaryCondition<Type,MeshType>>& bcs =
        this->mshLevel().mshField().boundaryConditions();

    labelVector lowerSlaveBoundary = zeroXYZ;
    labelVector upperSlaveBoundary = zeroXYZ;

    forAll(bcs, i)
    {
        const boundaryCondition<Type,MeshType>& bc = bcs[i];

        if
        (
            bc.boundaryOffsetDegree() == 1
         && bc.patch().type() == "boundary"
         && bc.slave()
        )
        {
            const label facei = faceNumber(bc.boundaryOffset());

            if (facei % 2 == 0)
            {
                lowerSlaveBoundary[facei/2] = 1;
            }
            else
            {
                upperSlaveBoundary[facei/2] = 1;
            }
        }
    }

    setActiveCells
    (
        briscola::cmptMax(fvMsh_.msh().lowerSlavePatch(), lowerSlaveBoundary),
        briscola::cmptMax(fvMsh_.msh().upperSlavePatch(), upperSlaveBoundary)
    );
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::allocate(const labelVector B)
{
    blockType::reAllocate(B);
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::transferData
(
    meshDirection<Type,MeshType>& D
)
{
    mshLevelPtr_ = D.mshLevelPtr_;

    l_ = D.l_;
    d_ = D.d_;

    blockType::transferData(D);

    D.mshLevelPtr_ = nullptr;
}

// Main constructor

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    meshLevel<Type,MeshType>& mshLevel,
    const fvMesh& fvMsh,
    const label l,
    const label d
)
:
    block<Type>(),
    fvMsh_(fvMsh),
    l_(l),
    d_(d),
    Ni_(fvMsh[l].N() + MeshType::padding[d]),
    mshLevelPtr_(&mshLevel)
{
    allocate(Ni_ + 2*unitXYZ);
    setActiveCells();
}

// Copy constructors

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const meshDirection<Type,MeshType>& D
)
:
    block<Type>(),
    fvMsh_(D.fvMsh_),
    l_(D.l_),
    d_(D.d_),
    S_(D.S_),
    E_(D.E_),
    N_(D.N_),
    Ni_(D.Ni_),
    mshLevelPtr_(nullptr)
{
    allocate(D.B().N());
    *this = D;
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const meshDirection<Type,MeshType>& D,
    const zero&
)
:
    block<Type>(),
    fvMsh_(D.fvMsh_),
    l_(D.l_),
    d_(D.d_),
    S_(D.S_),
    E_(D.E_),
    N_(D.N_),
    Ni_(D.Ni_),
    mshLevelPtr_(nullptr)
{
    allocate(D.B().N());
    *this = Zero;
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const meshDirection<Type,MeshType>& D,
    const Type& v
)
:
    block<Type>(),
    fvMsh_(D.fvMsh_),
    l_(D.l_),
    d_(D.d_),
    S_(D.S_),
    E_(D.E_),
    N_(D.N_),
    Ni_(D.Ni_),
    mshLevelPtr_(nullptr)
{
    allocate(D.B().N());
    *this = v;
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const tmp<meshDirection<Type,MeshType>>& tD
)
:
    block<Type>(),
    fvMsh_(tD->fvMsh_),
    l_(tD->l_),
    d_(tD->d_),
    S_(tD->S_),
    E_(tD->E_),
    N_(tD->N_),
    Ni_(tD->Ni_),
    mshLevelPtr_(nullptr)
{
    if (tD.isTmp())
    {
        meshDirection<Type,MeshType>& D =
            const_cast<meshDirection<Type,MeshType>&>(tD());

        transferData(D);
    }
    else
    {
        allocate(tD->B().N());
        *this = tD();
    }

    tD.clear();
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const tmp<meshDirection<Type,MeshType>>& tD,
    const zero&
)
:
    block<Type>(),
    fvMsh_(tD->fvMsh_),
    l_(tD->l_),
    d_(tD->d_),
    S_(tD->S_),
    E_(tD->E_),
    N_(tD->N_),
    Ni_(tD->Ni_),
    mshLevelPtr_(nullptr)
{
    if (tD.isTmp())
    {
        meshDirection<Type,MeshType>& D =
            const_cast<meshDirection<Type,MeshType>&>(tD());

        transferData(D);
    }
    else
    {
        allocate(tD->B().N());
    }

    *this = Zero;

    tD.clear();
}

template<class Type, class MeshType>
meshDirection<Type,MeshType>::meshDirection
(
    const tmp<meshDirection<Type,MeshType>>& tD,
    const Type& v
)
:
    block<Type>(),
    fvMsh_(tD->fvMsh_),
    l_(tD->l_),
    d_(tD->d_),
    S_(tD->S_),
    E_(tD->E_),
    N_(tD->N_),
    Ni_(tD->Ni_),
    mshLevelPtr_(nullptr)
{
    if (tD.isTmp())
    {
        meshDirection<Type,MeshType>& D =
            const_cast<meshDirection<Type,MeshType>&>(tD());

        transferData(D);
    }
    else
    {
        allocate(tD->B().N());
    }

    *this = v;

    tD.clear();
}

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
    Ni_(fvMsh[l].N() + MeshType::padding[d]),
    mshLevelPtr_(nullptr)
{
    allocate(Ni_ + 2*unitXYZ);
    setActiveCells();
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
    Ni_(fvMsh[l].N() + MeshType::padding[d]),
    mshLevelPtr_(nullptr)
{
    allocate(Ni_ + 2*unitXYZ);
    *this = Zero;
    setActiveCells();
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
    Ni_(fvMsh[l].N() + MeshType::padding[d]),
    mshLevelPtr_(nullptr)
{
    allocate(Ni_ + 2*unitXYZ);
    *this = v;
    setActiveCells();
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
void meshDirection<Type,MeshType>::operator=
(
    const tmp<meshDirection<Type,MeshType>>& tD
)
{
    if (tD.isTmp())
    {
        meshDirection<Type,MeshType>& D =
            const_cast<meshDirection<Type,MeshType>&>(tD());

        transferData(D);
    }
    else
    {
        *this = tD();
    }

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

}

}

}

#include "meshDirectionFunctions.C"
#include "meshDirectionStencilFunctions.C"
