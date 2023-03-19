#include "meshDirection.H"
#include "meshField.H"
#include "boundaryCondition.H"
#include "boundaryPartPatch.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::setInternalCells()
{
    const labelVector& padding = MeshType::padding[d_];

    const faceLabel slave = fvMsh_.msh().patchSlave();

    S_ = briscola::cmptMultiply(padding, slave.lower());

    E_ =
        this->B().N()
      - 2*unitXYZ
      - briscola::cmptMultiply(padding, slave.upper());

    N_ = E_ - S_;

    // By default, all internal cells are active cells

    Sa_ = S_;
    Ea_ = E_;
    Na_ = N_;
}

template<class Type, class MeshType>
void meshDirection<Type,MeshType>::updateActiveCells()
{
    const PtrList<boundaryCondition<Type,MeshType>>& bcs =
        this->mshLevel().mshField().boundaryConditions();

    faceLabel slaveBoundary(Zero);

    forAll(bcs, i)
    {
        const boundaryCondition<Type,MeshType>& bc = bcs[i];

        if
        (
            bc.boundaryOffsetDegree() == 1
         && bc.patch().type() == boundaryPartPatch::typeName
         && bc.slave()
        )
        {
            slaveBoundary[faceNumber(bc.boundaryOffset())] = 1;
        }
    }

    // If no boundary conditions are set the slave face label will revert to
    // that of the mesh, therewith setting Sa = S, Ea = E and Na = N

    const faceLabel slave =
        max(fvMsh_.msh().patchSlave(), slaveBoundary);

    const labelVector& padding = MeshType::padding[d_];

    Sa_ = briscola::cmptMultiply(padding, slave.lower());

    Ea_ =
        this->B().N()
      - 2*unitXYZ
      - briscola::cmptMultiply(padding, slave.upper());

    Na_ = Ea_ - Sa_;
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
    mshLevelPtr_(&mshLevel)
{
    allocate(fvMsh[l].N() + MeshType::padding[d] + 2*unitXYZ);

    setInternalCells();
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
    Sa_(D.Sa_),
    Ea_(D.Ea_),
    Na_(D.Na_),
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
    Sa_(D.Sa_),
    Ea_(D.Ea_),
    Na_(D.Na_),
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
    Sa_(D.Sa_),
    Ea_(D.Ea_),
    Na_(D.Na_),
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
    Sa_(tD->Sa_),
    Ea_(tD->Ea_),
    Na_(tD->Na_),
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
    Sa_(tD->Sa_),
    Ea_(tD->Ea_),
    Na_(tD->Na_),
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
    Sa_(tD->Sa_),
    Ea_(tD->Ea_),
    Na_(tD->Na_),
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
    mshLevelPtr_(nullptr)
{
    allocate(fvMsh[l].N() + MeshType::padding[d] + 2*unitXYZ);

    setInternalCells();
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
    mshLevelPtr_(nullptr)
{
    allocate(fvMsh[l].N() + MeshType::padding[d] + 2*unitXYZ);
    *this = Zero;

    setInternalCells();
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
    mshLevelPtr_(nullptr)
{
    allocate(fvMsh[l].N() + MeshType::padding[d] + 2*unitXYZ);
    *this = v;

    setInternalCells();
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
