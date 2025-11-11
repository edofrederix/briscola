#include "linearSystem.H"
#include "linearSystemAggregation.H"
#include "linearSystemFunctions.H"

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

    singular_ = sys.singular_;
    diagonal_ = sys.diagonal_;
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
    )
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
    setForcingMask();

    const label l = res.levelNum();
    const label d = res.directionNum();

    block<Type>& B = res.B();
    const labelVector shape = B.shape();

    // Restricted array pointers

    Type* const __restrict__ res_arr = B.begin();

    const Type* const __restrict__ x_arr = this->x()[l][d].B().begin();
    const Type* const __restrict__ b_arr = this->b()[l][d].B().begin();
    const label* const __restrict__ f_arr =
        this->forcingMask()[l][d].B().begin();

    // Reinterpret the matrix as a contiguous array of scalars

    const scalar* const __restrict__ A_arr =
        reinterpret_cast<const scalar*>(this->A()[l][d].B().begin());

    const faceLabel I = res.I();

    // Data in B is padded by ghosts

    const labelVector S = I.lower() + G*unitXYZ;
    const labelVector E = I.upper() + G*unitXYZ;
    const labelVector M = E - S;

    // Strides in i and j (data is contiguous in k)

    const label S_i = lin(S+unitX, shape) - lin(S, shape);
    const label S_j = lin(S+unitY, shape) - lin(S, shape);
    const label S_k = 1;

    // Jump after each line in k and plane in (j,k)

    const label J_k = lin(S+unitY, shape) - lin(S+unitZ*M.z(), shape);
    const label J_j = lin(S+unitX, shape) - lin(S+unitY*M.y(), shape);

    // Row product function

    linearSystemFun::rowProduct<SType,Type> P;

    // Compute residuals

    int c = lin(S, shape);

    if (NoMask || !this->x().immersedBoundaryConditions().size())
    {
        // Compute residual for each cell

        for (int i = 0; i < M.x(); i++)
        {
            for (int j = 0; j < M.y(); j++)
            {
                for (int k = 0; k < M.z(); k++)
                {
                    res_arr[c] =
                        b_arr[c]
                      - P(A_arr, x_arr, c, S_i, S_j, S_k);

                    // Jump to next cell
                    c++;
                }

                // Jump to next line
                c += J_k;
            }

            // Jump to next plane
            c += J_j;
        }
    }
    else
    {
        // Compute residual for each unmasked cell

        for (int i = 0; i < M.x(); i++)
        {
            for (int j = 0; j < M.y(); j++)
            {
                for (int k = 0; k < M.z(); k++)
                {
                    if (f_arr[c])
                    {
                        res_arr[c] = Zero;
                    }
                    else
                    {
                        res_arr[c] =
                            b_arr[c]
                          - P(A_arr, x_arr, c, S_i, S_j, S_k);
                    }

                    // Jump to next cell
                    c++;
                }

                // Jump to next line
                c += J_k;
            }

            // Jump to next plane
            c += J_j;
        }
    }
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
tmp<meshField<Type, MeshType>>
linearSystem<SType,Type,MeshType>::residual() const
{
    tmp<meshField<Type,MeshType>> tRes =
        meshField<Type,MeshType>::New
        (
            "residual",
            fvMsh_
        );

    this->residual<NoMask>(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::residual(const label l) const
{
    tmp<meshLevel<Type,MeshType>> tRes =
        meshLevel<Type,MeshType>::New(fvMsh_,l);

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
    tmp<meshDirection<Type,MeshType>> tRes =
        meshDirection<Type,MeshType>::New(fvMsh_,l,d);

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
    setForcingMask();

    const label l = eval.levelNum();
    const label d = eval.directionNum();

    block<Type>& B = eval.B();
    const labelVector shape = B.shape();

    // Restricted array pointers

    Type* const __restrict__ eval_arr = B.begin();

    const Type* const __restrict__ x_arr = this->x()[l][d].B().begin();
    const Type* const __restrict__ b_arr = this->b()[l][d].B().begin();
    const label* const __restrict__ f_arr =
        this->forcingMask()[l][d].B().begin();

    const scalar* const __restrict__ icv_arr =
        this->fvMsh_.template
        metrics<MeshType>().inverseCellVolumes()[l][d].B().begin();

    // Reinterpret the matrix as a contiguous array of scalars

    const scalar* const __restrict__ A_arr =
        reinterpret_cast<const scalar*>(this->A()[l][d].B().begin());

    const faceLabel I = eval.I();

    // Data in B is padded by ghosts

    const labelVector S = I.lower() + G*unitXYZ;
    const labelVector E = I.upper() + G*unitXYZ;
    const labelVector M = E - S;

    // Strides in i and j (data is contiguous in k)

    const label S_i = lin(S+unitX, shape) - lin(S, shape);
    const label S_j = lin(S+unitY, shape) - lin(S, shape);
    const label S_k = 1;

    // Jump after each line in k and plane in (j,k)

    const label J_k = lin(S+unitY, shape) - lin(S+unitZ*M.z(), shape);
    const label J_j = lin(S+unitX, shape) - lin(S+unitY*M.y(), shape);

    // Row product function

    linearSystemFun::rowProduct<SType,Type> P;

    // Compute evaluations

    int c = lin(S, shape);

    if (NoMask || !this->x().immersedBoundaryConditions().size())
    {
        // Compute evaluation for each cell

        for (int i = 0; i < M.x(); i++)
        {
            for (int j = 0; j < M.y(); j++)
            {
                for (int k = 0; k < M.z(); k++)
                {
                    eval_arr[c] =
                        (
                            P(A_arr, x_arr, c, S_i, S_j, S_k)
                          - b_arr[c]
                        )
                      * icv_arr[c];

                    // Jump to next cell
                    c++;
                }

                // Jump to next line
                c += J_k;
            }

            // Jump to next plane
            c += J_j;
        }
    }
    else
    {
        // Compute evaluation for each unmasked cell

        for (int i = 0; i < M.x(); i++)
        {
            for (int j = 0; j < M.y(); j++)
            {
                for (int k = 0; k < M.z(); k++)
                {
                    if (f_arr[c])
                    {
                        eval_arr[c] = Zero;
                    }
                    else
                    {
                        eval_arr[c] =
                            (
                                P(A_arr, x_arr, c, S_i, S_j, S_k)
                              - b_arr[c]
                            )
                          * icv_arr[c];
                    }

                    // Jump to next cell
                    c++;
                }

                // Jump to next line
                c += J_k;
            }

            // Jump to next plane
            c += J_j;
        }
    }
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
tmp<meshField<Type, MeshType>>
linearSystem<SType,Type,MeshType>::evaluate() const
{
    tmp<meshField<Type,MeshType>> tEval =
        meshField<Type,MeshType>::New
        (
            "evaluate",
            fvMsh_
        );

    this->evaluate<NoMask>(tEval.ref());

    return tEval;
}

template<class SType, class Type, class MeshType>
template<bool NoMask>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::evaluate(const label l) const
{
    tmp<meshLevel<Type,MeshType>> tEval =
        meshLevel<Type,MeshType>::New(fvMsh_,l);

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
    tmp<meshDirection<Type,MeshType>> tRes =
        meshDirection<Type,MeshType>::New(fvMsh_,l,d);

    this->evaluate<NoMask>(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::eliminateGhosts()
{
    if (!xPtr_->boundaryConditions().size())
        xPtr_->addBoundaryConditions();

    forAll(A_, l)
        forAll(xPtr_->boundaryConditions(), i)
            xPtr_->boundaryConditions()[i].eliminateGhosts(*this, l);
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

        #ifdef NO_BLOCK_ZERO_INIT
        f = Zero;
        #endif

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

    this->singular_ = sys.singular();
    this->diagonal_ = sys.diagonal();
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

    this->singular_.clear();
    this->diagonal_.clear();
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

    this->singular_.clear();
    this->diagonal_.clear();
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

    this->singular_.clear();
    this->diagonal_.clear();
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

    this->singular_.clear();
    this->diagonal_.clear();
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

    singular_.clear();
    diagonal_.clear();
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

    this->singular_.clear();
    this->diagonal_.clear();
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

    this->singular_.clear();
    this->diagonal_.clear();
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

    this->singular_.clear();
    this->diagonal_.clear();
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
    restrict(field);

    this->A() *= field;
    this->b() *= field;

    this->singular_.clear();
    this->diagonal_.clear();

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
    restrict(field);

    this->A() /= field;
    this->b() /= field;

    this->singular_.clear();
    this->diagonal_.clear();

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
