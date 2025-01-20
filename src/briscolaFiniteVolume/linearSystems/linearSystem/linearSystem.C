#include "linearSystem.H"
#include "linearSystemAggregation.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::transfer
(
    linearSystem<SType,Type,MeshType>& sys
)
{
    A_.transfer(sys.A_);
    b_.transfer(sys.b_);
}

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
            x.fvMsh().time().timeName(),
            x.fvMsh().time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    refCount(),
    fvMsh_(x.fvMsh()),
    x_(x),
    A_
    (
        IOobject::groupName("A", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
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
            sys.fvMsh_.time().timeName(),
            sys.fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    refCount(),
    fvMsh_(sys.fvMsh_),
    x_(sys.x_),
    A_
    (
        IOobject::groupName("A", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
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
            tSys->fvMsh_.time().timeName(),
            tSys->fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    refCount(),
    fvMsh_(tSys->fvMsh_),
    x_(tSys->x_),
    A_
    (
        IOobject::groupName("A", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
{
    if (tSys.isTmp())
    {
        linearSystem<SType,Type,MeshType>& sys =
            const_cast<linearSystem<SType,Type,MeshType>&>(tSys());

        transfer(sys);
    }
    else
    {
        *this = tSys();
    }

    if (tSys.isTmp())
        tSys.clear();
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
            sys.fvMsh_.time().timeName(),
            sys.fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    refCount(),
    fvMsh_(sys.fvMsh_),
    x_(x),
    A_
    (
        IOobject::groupName("A", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
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
            tSys->fvMsh_.time().timeName(),
            tSys->fvMsh_.time(),
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            registerObject
        )
    ),
    refCount(),
    fvMsh_(tSys->fvMsh_),
    x_(x),
    A_
    (
        IOobject::groupName("A", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE,
        false,
        true
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    )
{
    if (tSys.isTmp())
    {
        linearSystem<SType,Type,MeshType>& sys =
            const_cast<linearSystem<SType,Type,MeshType>&>(tSys());

        transfer(sys);
    }
    else
    {
        *this = tSys();
    }

    if (tSys.isTmp())
        tSys.clear();
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::~linearSystem()
{}

template<class SType, class Type, class MeshType>
List<bool> linearSystem<SType,Type,MeshType>::singular(const bool clearCache)
{
    if (singular_.size() == 0 || clearCache)
    {
        singular_.setSize(MeshType::numberOfDirections);
        singular_ = true;

        forAll(singular_, d)
        {
            // We prefer a level with at least three cells in one direction to
            // avoid looking at only constrained cells

            label l = A_.size() - 1;
            while (l > 0 && cmptMax(A_[l][d].N()) < 3) l--;

            const meshDirection<SType,MeshType>& A = A_[l][d];

            forAllCells(A, i, j, k)
                if (Foam::mag(rowSum(A,i,j,k)) > 1e-8)
                    goto nope;

            continue;
            nope: singular_[d] = false;
        }

        reduce(singular_, andOp<List<bool>>());
    }

    return singular_;
}

template<class SType, class Type, class MeshType>
List<bool> linearSystem<SType,Type,MeshType>::diagonal(const bool clearCache)
{
    if (diagonal_.size() == 0 || clearCache)
    {
        diagonal_.setSize(MeshType::numberOfDirections);
        diagonal_ = true;

        forAll(diagonal_, d)
        {
            // We prefer a level with at least three cells in one direction to
            // avoid looking at only constrained cells

            label l = A_.size() - 1;
            while (l > 0 && cmptMax(A_[l][d].N()) < 3) l--;

            const meshDirection<SType,MeshType>& A = A_[l][d];

            forAllCells(A, i, j, k)
                if (A(i,j,k) != SType(diagStencil(A(i,j,k).center())))
                    goto nope;

            continue;
            nope: diagonal_[d] = false;
        }

        reduce(diagonal_, andOp<List<bool>>());
    }

    return diagonal_;
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
void linearSystem<SType,Type,MeshType>::residual
(
    meshField<Type, MeshType>& res
) const
{
    forAll(res, l)
        this->residual<NoMask>(res[l]);
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
void linearSystem<SType,Type,MeshType>::residual
(
    meshLevel<Type,MeshType>& res
) const
{
    forAll(res, d)
        this->residual<NoMask>(res[d]);
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
void linearSystem<SType,Type,MeshType>::residual
(
    meshDirection<Type,MeshType>& res
) const
{
    const label l = res.levelNum();
    const label d = res.directionNum();

    const meshDirection<SType,MeshType>& A = this->A()[l][d];
    const meshDirection<Type,MeshType>& x = this->x()[l][d];
    const meshDirection<Type,MeshType>& b = this->b()[l][d];

    if (diagonal_.size() && diagonal_[d])
    {
        forAllCells(res, i, j, k)
            res(i,j,k) = b(i,j,k) - A(i,j,k).center()*x(i,j,k);
    }
    else
    {
        forAllCells(res, i, j, k)
            res(i,j,k) = b(i,j,k) - rowProduct(A,x,i,j,k);
    }

    if (!NoMask && x_.immersedBoundaryConditions().size())
    {
        setForcingMask();

        const meshDirection<label,MeshType>& f = forcingMask_()[l][d];

        forAllCells(res,i,j,k)
            if (f(i,j,k))
                res(i,j,k) = Zero;
    }
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
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

    this->residual<NoMask>(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::residual(const label l) const
{
    tmp<meshLevel<Type,MeshType>> tRes
    (
        new meshLevel<Type,MeshType>(fvMsh_,l)
    );

    this->residual<NoMask>(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
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

    this->residual<NoMask>(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
void linearSystem<SType,Type,MeshType>::evaluate
(
    meshField<Type, MeshType>& eval
) const
{
    forAll(eval, l)
        this->evaluate<NoMask>(eval[l]);
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
void linearSystem<SType,Type,MeshType>::evaluate
(
    meshLevel<Type,MeshType>& eval
) const
{
    forAll(eval, d)
        this->evaluate<NoMask>(eval[d]);
}


template<class SType, class Type, class MeshType>
template<bool NoMask>
void linearSystem<SType,Type,MeshType>::evaluate
(
    meshDirection<Type,MeshType>& eval
) const
{
    const label l = eval.levelNum();
    const label d = eval.directionNum();

    const meshDirection<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes()[l][d];

    const meshDirection<SType,MeshType>& A = this->A()[l][d];
    const meshDirection<Type,MeshType>& x = this->x()[l][d];
    const meshDirection<Type,MeshType>& b = this->b()[l][d];

    if (diagonal_.size() && diagonal_[d])
    {
        forAllCells(eval, i, j, k)
            eval(i,j,k) = (A(i,j,k).center()*x(i,j,k) - b(i,j,k))/cv(i,j,k);
    }
    else
    {
        forAllCells(eval, i, j, k)
            eval(i,j,k) = (rowProduct(A,x,i,j,k) - b(i,j,k))/cv(i,j,k);
    }

    if (!NoMask && x_.immersedBoundaryConditions().size())
    {
        setForcingMask();

        const meshDirection<label,MeshType>& f = forcingMask_()[l][d];

        forAllCells(eval,i,j,k)
            if (f(i,j,k))
                eval(i,j,k) = Zero;
    }
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
tmp<meshField<Type, MeshType>>
linearSystem<SType,Type,MeshType>::evaluate() const
{
    tmp<meshField<Type,MeshType>> tEval
    (
        new meshField<Type,MeshType>
        (
            "evaluate",
            fvMsh_
        )
    );

    this->evaluate<NoMask>(tEval.ref());

    return tEval;
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::evaluate(const label l) const
{
    tmp<meshLevel<Type,MeshType>> tEval
    (
        new meshLevel<Type,MeshType>(fvMsh_,l)
    );

    this->evaluate<NoMask>(tEval.ref());

    return tEval;
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
tmp<meshDirection<Type,MeshType>>
linearSystem<SType,Type,MeshType>::evaluate
(
    const label l,
    const label d
) const
{
    tmp<meshDirection<Type,MeshType>> tRes
    (
        new meshDirection<Type,MeshType>(fvMsh_,l,d)
    );

    this->evaluate<NoMask>(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::eliminateGhosts()
{
    if (!x_.boundaryConditions().size())
        x_.addBoundaryConditions();

    forAll(A_, l)
        forAll(x_.boundaryConditions(), i)
            x_.boundaryConditions()[i].eliminateGhosts(*this, l);
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::setForcingMask()
{
    if (forcingMask_.empty())
    {
        forcingMask_.set
        (
            new meshField<label,MeshType>
            (
                IOobject::groupName(x_.name(), "forcingMask"),
                fvMsh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                true,
                true
            )
        );

        meshField<label,MeshType>& f = forcingMask_();

        f = Zero;

        forAll(x_.immersedBoundaryConditions(), i)
            if (x_.immersedBoundaryConditions()[i].forcingMaskPtr())
                f += x_.immersedBoundaryConditions()[i].forcingMask();

        f = min(f,1);
        f.correctBoundaryConditions();
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

        transfer(sys);
    }
    else
    {
        *this = tSys();
    }

    if (tSys.isTmp())
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
    const Type& v
)
{
    const meshField<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes();

    this->A() = diagStencil(0.0);
    this->b() = -cv*v;
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
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator=
(
    const linearSystem<SType2,Type,MeshType>& sys
)
{
    this->A() = sys.A();
    this->b() = sys.b();
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
    const SType2& v
)
{
    this->b() -= v;
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
    const SType2& v
)
{
    this->b() += v;
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
    const_cast<meshField<scalar,MeshType>&>(field).restrict();

    this->A() *= field;
    this->b() *= field;

    const_cast<meshField<scalar,MeshType>&>(field).makeShallow();
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
    const_cast<meshField<scalar,MeshType>&>(field).restrict();

    this->A() /= field;
    this->b() /= field;

    const_cast<meshField<scalar,MeshType>&>(field).makeShallow();
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
