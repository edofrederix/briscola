#include "linearSystem.H"
#include "Row.H"

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
            // We need at least three cells in one direction to avoid looking at
            // only constrained cells

            label l = A_.size() - 1;
            while (cmptMax(A_[l][d].N()) < 3) l--;

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
            // We need at least three cells in one direction to avoid looking at
            // only constrained cells

            label l = A_.size() - 1;
            while (cmptMax(A_[l][d].N()) < 3) l--;

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
void linearSystem<SType,Type,MeshType>::residual
(
    meshField<Type, MeshType>& res
)
{
    forAll(res, l)
        this->residual(res[l]);
}

template<class SType, class Type, class MeshType>
tmp<meshField<Type, MeshType>>
linearSystem<SType,Type,MeshType>::residual()
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
)
{
    forAll(res, d)
        this->residual(res[d]);
}

template<class SType, class Type, class MeshType>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::residual(const label l)
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
)
{
    const label l = res.levelNum();
    const label d = res.directionNum();

    const meshDirection<SType,MeshType>& A = this->A()[l][d];
    const meshDirection<Type,MeshType>& x = this->x()[l][d];
    const meshDirection<Type,MeshType>& b = this->b()[l][d];

    res = Zero;

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

    setIBMForcingMask();

    if (!x_.immersedBoundaryConditions().empty())
    {
        forAllCells(res,i,j,k)
        {
            if (IBMForcingMask_()(l,d,i,j,k))
            {
                res(i,j,k) = Zero;
            }
        }
    }
}

template<class SType, class Type, class MeshType>
tmp<meshDirection<Type,MeshType>>
linearSystem<SType,Type,MeshType>::residual
(
    const label l,
    const label d
)
{
    tmp<meshDirection<Type,MeshType>> tRes
    (
        new meshDirection<Type,MeshType>(fvMsh_,l,d)
    );

    this->residual(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::evaluate
(
    meshField<Type, MeshType>& eval
) const
{
    forAll(eval, l)
        this->evaluate(eval[l]);
}

template<class SType, class Type, class MeshType>
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

    this->evaluate(tEval.ref());

    return tEval;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::evaluate
(
    meshLevel<Type,MeshType>& eval
) const
{
    forAll(eval, d)
        this->evaluate(eval[d]);
}

template<class SType, class Type, class MeshType>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::evaluate(const label l) const
{
    tmp<meshLevel<Type,MeshType>> tEval
    (
        new meshLevel<Type,MeshType>(fvMsh_,l)
    );

    this->evaluate(tEval.ref());

    return tEval;
}

template<class SType, class Type, class MeshType>
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

    eval = Zero;

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
}

template<class SType, class Type, class MeshType>
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

    this->evaluate(tRes.ref());

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
void linearSystem<SType,Type,MeshType>::setIBMForcingMask()
{
    if (IBMForcingMask_.empty())
    {
        IBMForcingMask_.set
        (
            new meshField<label,MeshType>
            (
                IOobject::groupName(x_.name(),"IBMForcingMask"),
                fvMsh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false,
                true
            )
        );

        IBMForcingMask_() = Zero;

        if (!x_.immersedBoundaryConditions().empty())
        {
            forAll(x_.immersedBoundaryConditions(), ib)
            {
                forAllCells(IBMForcingMask_(),l,d,i,j,k)
                {
                    if
                    (
                        x_.immersedBoundaryConditions()[ib]
                            .forcingPoints()(l,d,i,j,k)
                    )
                    {
                        IBMForcingMask_()(l,d,i,j,k) = 1;
                    }
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
    this->A() *= field;
    this->b() *= field;
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
    this->A() /= field;
    this->b() /= field;
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
void writeToFile
(
    linearSystem<SType,Type,MeshType>& sys,
    const fileName file,
    const label l = 0
)
{
    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        const meshDirection<SType,MeshType>& A = sys.A()[l][d];
        const meshDirection<Type,MeshType>& b = sys.b()[l][d];

        const labelVector N = A.N();

        List<labelVector> shapes(Pstream::nProcs());
        shapes[Pstream::myProcNo()] = N;

        Pstream::gatherList(shapes);

        List<List<Row<stencil,Type>>> data
        (
            Pstream::nProcs()
        );

        data[Pstream::myProcNo()].setSize(cmptProduct(N));

        int l = 0;
        forAllCells(A, i, j, k)
            data[Pstream::myProcNo()][l++] =
                Row<stencil,Type>
                (
                    fullStencil<MeshType>(A,i,j,k),
                    b(i,j,k)
                );

        Pstream::gatherList(data);

        if (Pstream::master())
        {
            OFstream F
            (
                file
              + fileName
                (
                    MeshType::numberOfDirections > 1
                  ? ("_" + Foam::name(d)) : ""
                )
            );

            label n = 0;
            forAll(data, i)
                n += data[i].size();

            labelList offsets(Pstream::nProcs());

            offsets[0] = 0;
            for (int proc = 1; proc < Pstream::nProcs(); proc++)
                offsets[proc] = offsets[proc-1] + cmptProduct(shapes[proc-1]);

            const decomposition& decomp = sys.fvMsh().msh().decomp();

            forAll(data, proc)
            {
                const labelVector N = shapes[proc];

                forAll(data[proc], i)
                {
                    const stencil Ai(data[proc][i].stencil());
                    const Type bi(data[proc][i].source());

                    const labelVector ijk
                    (
                        (i/N.y()/N.z()) % N.x(),
                        (i/N.z()) % N.y(),
                        i % N.z()
                    );

                    faceLabel neigh(-faceLabel::one);

                    for (int f = 0; f < 6; f++)
                    if (Ai[f+1] != 0.0)
                    {
                        labelVector ijkn(ijk + faceOffsets[f]);

                        if
                        (
                            briscola::cmptMax
                            (
                                briscola::cmptMin(ijkn, N-unitXYZ),
                                zeroXYZ
                            )
                          == ijkn
                        )
                        {
                            // Local cell

                            neigh[f] =
                                offsets[proc]
                              + ijkn.x()*N.y()*N.z()
                              + ijkn.y()*N.z()
                              + ijkn.z();
                        }
                        else
                        {
                            // Neighbor cell

                            const label neighProc =
                                decomp.faceNeighborsPerProc()[proc][f];

                            if (neighProc > -1)
                            {
                                const labelTensor T =
                                    decomp.faceTsPerProc()[proc][f];

                                const labelVector Nn = shapes[neighProc];
                                const labelVector TNn = cmptMag(T & Nn);

                                ijkn[f/2] = ijkn[f/2] < 0 ? TNn[f/2]-1 : 0;

                                ijkn = (T.T() & ijkn) + Nn;

                                for (int j = 0; j < 3; j++)
                                    ijkn[j] = (ijkn[j] % Nn[j]);

                                neigh[f] =
                                    offsets[neighProc]
                                  + ijkn.x()*Nn.y()*Nn.z()
                                  + ijkn.y()*Nn.z()
                                  + ijkn.z();
                            }
                        }
                    }

                    // Write

                    for (int j = 0; j < n; j++)
                    {
                        if (offsets[proc]+i == j)
                        {
                            F << Ai[0] << " ";
                        }
                        else
                        {
                            bool found = false;
                            for (int f = 0; f < 6; f++)
                            {
                                if (neigh[f] == j)
                                {
                                    F << Ai[f+1] << " ";
                                    found = true;
                                    break;
                                }
                            }

                            if (!found)
                                F << 0.0 << " ";
                        }
                    }

                    F << bi << nl;
                }
            }
        }
    }
}

}

}

}
