#include "linearSystem.H"
#include "linearSystemAggregation.H"
#include "linearSystemFunctions.H"

#include "boundaryConditionSelector.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

// Private functions

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::transfer
(
    linearSystem<SType,Type,MeshType>& sys
)
{
    A_.transfer(sys.A_);
    b_.transfer(sys.b_);

    symmetric_ = sys.symmetric_;
    diagonal_ = sys.diagonal_;
    eliminated_ = sys.eliminated_;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::correctEliminatedBoundaryMasks()
{
    eliminatedBoundaryMasks_.resize(fvMsh_.size());
    eliminatedBoundaryMasks_ = Zero;

    const meshField<Type,MeshType>& x = *xPtr_;

    if (eliminated_)
        forAll(x, l)
            forAll(x[l].boundaryConditions(), i)
                if (x[l].boundaryConditions()[i].offsetDegree() == 1)
                    eliminatedBoundaryMasks_[l]
                        [faceNumber(x[l].boundaryConditions()[i].offset())] =
                            x[l].boundaryConditions()[i].eliminated();
}

// Constructors

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    const word name,
    meshField<Type, MeshType>& x,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            name,
            x.fvMsh().time().name(),
            x.fvMsh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(x.fvMsh()),
    xPtr_(&x),
    A_
    (
        IOobject::groupName(IOobject::name(), "A"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName(IOobject::name(), "b"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    symmetric_(true),
    diagonal_(true),
    eliminated_(false)
{}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    const IOobject& io,
    const fvMesh& fvMsh
)
:
    regIOobject(io),
    cachedRefCount(),
    fvMsh_(fvMsh),
    xPtr_(nullptr),
    A_
    (
        IOobject::groupName(IOobject::name(), "A"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName(IOobject::name(), "b"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    symmetric_(true),
    diagonal_(true),
    eliminated_(false)
{}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    const linearSystem<SType,Type,MeshType>& sys,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            sys.name(),
            sys.fvMsh_.time().name(),
            sys.fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(sys.fvMsh_),
    xPtr_(sys.xPtr_),
    A_
    (
        IOobject::groupName(IOobject::name(), "A"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName(IOobject::name(), "b"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
{
    *this = sys;
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            tSys->name(),
            tSys->fvMsh_.time().name(),
            tSys->fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(tSys->fvMsh_),
    xPtr_(tSys->xPtr_),
    A_
    (
        IOobject::groupName(IOobject::name(), "A"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName(IOobject::name(), "b"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
{
    if (tSys.isTmp() && tSys->unique())
    {
        linearSystem<SType,Type,MeshType>& sys =
            const_cast<linearSystem<SType,Type,MeshType>&>(tSys());

        transfer(sys);

        tSys.clear();
    }
    else
    {
        *this = tSys();
    }
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    const linearSystem<SType,Type,MeshType>& sys,
    meshField<Type, MeshType>& x,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            sys.name(),
            sys.fvMsh_.time().name(),
            sys.fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(sys.fvMsh_),
    xPtr_(&x),
    A_
    (
        IOobject::groupName(IOobject::name(), "A"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName(IOobject::name(), "b"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
{
    *this = sys;
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::linearSystem
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys,
    meshField<Type, MeshType>& x,
    const bool registerObject
)
:
    regIOobject
    (
        IOobject
        (
            tSys->name(),
            tSys->fvMsh_.time().name(),
            tSys->fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    cachedRefCount(),
    fvMsh_(tSys->fvMsh_),
    xPtr_(&x),
    A_
    (
        IOobject::groupName(IOobject::name(), "A"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName(IOobject::name(), "b"),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
{
    if (tSys.isTmp() && tSys->unique())
    {
        linearSystem<SType,Type,MeshType>& sys =
            const_cast<linearSystem<SType,Type,MeshType>&>(tSys());

        transfer(sys);

        tSys.clear();
    }
    else
    {
        *this = tSys();
    }
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::~linearSystem()
{}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::eliminateGhosts()
{
    meshField<Type,MeshType>& x = *xPtr_;

    x.addBoundaryConditions();

    // Loop over the levels of x, so that if x is shallow, only the top level is
    // handled.

    forAll(x, l)
        forAll(x[l].boundaryConditions(), i)
            x[l].boundaryConditions()[i].eliminateGhosts(*this);

    eliminated_ = true;
    correctEliminatedBoundaryMasks();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::setForcingMask()
{
    if (forcingMask_.empty())
    {
        forcingMask_.set
        (
            meshField<label,MeshType>::New
            (
                IOobject::groupName(xPtr_->name(), "forcingMask"),
                fvMsh_
            ).ptr()
        );

        forcingMask_->makeDeep();

        meshField<label,MeshType>& f = forcingMask_();

        f = Zero;

        if (xPtr_->immersedBoundaryConditions().size())
        {
            forAll(xPtr_->immersedBoundaryConditions(), i)
                if (xPtr_->immersedBoundaryConditions()[i].forcingMaskPtr())
                    f += xPtr_->immersedBoundaryConditions()[i].forcingMask();

            f = min(f,1);
            f.correctBoundaryConditions();
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

    this->symmetric_ = sys.symmetric();
    this->diagonal_ = sys.diagonal();
    this->eliminated_ = sys.eliminated();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator=
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    if (tSys.isTmp() && tSys->unique())
    {
        linearSystem<SType,Type,MeshType>& sys =
            const_cast<linearSystem<SType,Type,MeshType>&>(tSys());

        transfer(sys);

        tSys.clear();
    }
    else
    {
        *this = tSys();
    }
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const linearSystem<SType,Type,MeshType>& sys
)
{
    this->A() += sys.A();
    this->b() += sys.b();

    symmetric_ = symmetric_ && sys.symmetric();
    diagonal_ = diagonal_ && sys.diagonal();
    eliminated_ = eliminated_ && sys.eliminated();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    *this += tSys();

    if (tSys.isTmp())
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

    symmetric_ = symmetric_ && sys.symmetric();
    diagonal_ = diagonal_ && sys.diagonal();
    eliminated_ = eliminated_ && sys.eliminated();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    *this -= tSys();

    if (tSys.isTmp())
        tSys.clear();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator=
(
    const zero
)
{
    this->A() = Zero;
    this->b() = Zero;

    symmetric_ = true;
    diagonal_ = true;
    eliminated_ = false;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator=
(
    const Type& v
)
{
    const meshField<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes();

    this->A() = Zero;
    this->b() = -cv*v;

    symmetric_ = true;
    diagonal_ = true;
    eliminated_ = false;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator=
(
    const meshField<Type,MeshType>& field
)
{
    const meshField<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes();

    this->A() = Zero;
    this->b() = -cv*field;

    symmetric_ = true;
    diagonal_ = true;
    eliminated_ = false;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator=
(
    const tmp<meshField<Type,MeshType>>& tField
)
{
    *this = tField();

    if (tField.isTmp())
        tField.clear();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const Type& v
)
{
    const meshField<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes();

    this->b() -= cv*v;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const List<Type>& v
)
{
    const meshField<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes();

    this->b() -= cv*v;
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

    if (tField.isTmp())
        tField.clear();
}

template<class SType, class Type, class MeshType>
template<class Type2>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const Type2& v
)
{
    this->b() -= v;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const Type& v
)
{
    const meshField<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes();

    this->b() += cv*v;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const List<Type>& v
)
{
    const meshField<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes();

    this->b() += cv*v;
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

    if (tField.isTmp())
        tField.clear();
}

template<class SType, class Type, class MeshType>
template<class Type2>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const Type2& v
)
{
    this->b() += v;
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

    symmetric_ = sys.symmetric();
    diagonal_ = sys.diagonal();
    eliminated_ = sys.eliminated();
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator=
(
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys
)
{
    *this = tSys();

    if (tSys.isTmp())
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

    symmetric_ = symmetric_ && sys.symmetric();
    diagonal_ = diagonal_ && sys.diagonal();
    eliminated_ = eliminated_ && sys.eliminated();
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys
)
{
    *this += tSys();

    if (tSys.isTmp())
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

    symmetric_ = symmetric_ && sys.symmetric();
    diagonal_ = diagonal_ && sys.diagonal();
    eliminated_ = eliminated_ && sys.eliminated();
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys
)
{
    *this -= tSys();

    if (tSys.isTmp())
        tSys.clear();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator*=
(
    const scalar s
)
{
    this->A() *= s;
    this->b() *= s;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator*=
(
    const scalarList& s
)
{
    this->A() *= s;
    this->b() *= s;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator*=
(
    const meshField<scalar,MeshType>& field
)
{
    const bool shallow = field.shallow();

    if (shallow)
        restrict(field);

    this->A() *= field;
    this->b() *= field;

    if (shallow)
        collapse(field);
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator*=
(
    const tmp<meshField<scalar,MeshType>>& tField
)
{
    *this *= tField();

    if (tField.isTmp())
        tField.clear();
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator/=
(
    const scalar s
)
{
    this->A() /= s;
    this->b() /= s;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator/=
(
    const scalarList& s
)
{
    this->A() /= s;
    this->b() /= s;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator/=
(
    const meshField<scalar,MeshType>& field
)
{
    const bool shallow = field.shallow();

    if (shallow)
        restrict(field);

    this->A() /= field;
    this->b() /= field;

    if (shallow)
        collapse(field);
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator/=
(
    const tmp<meshField<scalar,MeshType>>& tField
)
{
    *this /= tField();

    if (tField.isTmp())
        tField.clear();
}

template<class SType, class Type, class MeshType>
bool linearSystem<SType,Type,MeshType>::writeLevel(Ostream& os, const label l)
const
{
    // Aggregated linear system at master
    const linearSystemAggregation<SType,Type,MeshType> lsa(*this,l,1);

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        List<List<SType>> rows;
        lsa.rowCoeffs(rows, *this, d);

        List<Type> rhs;
        lsa.rhsSource(rhs, *this, d);

        if (Pstream::master())
        {
            const label n = rhs.size();

            os << n << " " << n << nl;

            const auto& colNums = lsa.colNums()[d];

            int c = 0;
            forAll(rows, proci)
            {
                forAll(rows[proci], row)
                {
                    const auto& cols = colNums[proci][row];

                    for (int col = 0; col < n; col++)
                    {
                        const label i = findIndex(cols, col);

                        if (i > -1 && rows[proci][row][i] != 0)
                        {
                            os << rows[proci][row][i] << " ";
                        }
                        else
                        {
                            os << "0 ";
                        }
                    }

                    for (int i = 0; i < pTraits<Type>::nComponents; i++)
                        os << scalar_cast(&rhs[c])[i] << " ";

                    os << nl;

                    c++;
                }
            }
        }
    }

    return os;
}

}

}

}
