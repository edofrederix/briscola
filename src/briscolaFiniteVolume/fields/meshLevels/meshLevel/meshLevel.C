#include "meshLevel.H"
#include "meshField.H"
#include "fvMesh.H"

#include "boundaryCondition.H"
#include "boundaryConditionSelector.H"
#include "immersedBoundaryCondition.H"

#include "patchBoundary.H"
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
            new meshDirection<Type,MeshType>(*this, d)
        );
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::setLevelPointers()
{
    forAll(*this, d)
        listType::operator[](d).levelPtr_ = this;
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::transfer
(
    meshLevel<Type,MeshType>& L
)
{
    listType::transfer(L);

    fieldPtr_ = L.fieldPtr_;
    L.fieldPtr_ = nullptr;

    setLevelPointers();
}

// Main constructor

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    meshField<Type,MeshType>& field,
    const label l
)
:
    FastPtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(field.fvMsh()),
    l_(l),
    fieldPtr_(&field)
{
    allocate();
}

// Copy constructors

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel(const meshLevel<Type,MeshType>& L)
:
    FastPtrList<meshDirection<Type,MeshType>>(L),
    refCount(),
    fvMsh_(L.fvMsh_),
    l_(L.levelNum()),
    fieldPtr_(nullptr)
{
    setLevelPointers();
    transferBoundaryConditions(L);
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel(const tmp<meshLevel<Type,MeshType>>& tL)
:
    FastPtrList<meshDirection<Type,MeshType>>
    (
        const_cast<meshLevel<Type,MeshType>&>(tL()),
        tL.isTmp()
    ),
    refCount(),
    fvMsh_(tL->fvMsh_),
    l_(tL->levelNum()),
    fieldPtr_(nullptr)
{
    setLevelPointers();
    transferBoundaryConditions(tL());

    if (tL.isTmp())
        tL.clear();
}

// Standalone constructors

template<class Type, class MeshType>
meshLevel<Type,MeshType>::meshLevel
(
    const fvMesh& fvMsh,
    const label l
)
:
    FastPtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(fvMsh),
    l_(l),
    fieldPtr_(nullptr)
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
    FastPtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(fvMsh),
    l_(l),
    fieldPtr_(nullptr)
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
    FastPtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(fvMsh),
    l_(l),
    fieldPtr_(nullptr)
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
    FastPtrList<meshDirection<Type,MeshType>>(),
    refCount(),
    fvMsh_(fvMsh),
    l_(l),
    fieldPtr_(nullptr)
{
    allocate();
    *this = v;
}

template<class Type, class MeshType>
meshLevel<Type,MeshType>::~meshLevel()
{}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::addBoundaryConditions()
{
    if (boundaryConditions_.empty())
    {
        boundaryConditions_.setSize(lvl().boundaries().size());

        singular_ = true;

        forAll(lvl().boundaries(), i)
        {
            boundaryConditions_.set
            (
                i,
                boundaryCondition<Type,MeshType>::New
                (
                    *this,
                    lvl().boundaries()[i]
                )
            );

            const boundaryConditionBaseType bcType =
                boundaryConditions_[i].baseType();

            if (bcType == DIRICHLETBC || bcType == ROBINBC)
                singular_ = false;
        }

        reduce(singular_, andOp<bool>(), Pstream::msgType(), lvl().comms());
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::transferBoundaryConditions
(
    const meshLevel<Type,MeshType>& L
)
{
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

    PtrList<boundaryCondition<Type,MeshType>> list
    (
        L.boundaryConditions(),
        *this
    );

    boundaryConditions_.clear();
    boundaryConditions_.transfer(list);
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctBoundaryConditions()
{
    if (fieldPtr_)
    {
        addBoundaryConditions();

        // Immersed boundary conditions must be corrected first, so that cell
        // values near the boundaries can be copied to neighbors. By default,
        // only immersed boundary conditions on the first level are corrected.
        // If the immersed boundary conditions must be corrected on another
        // level, then the correctImmersedBoundaryConditions() must be called
        // explicitly on that level.

        if (l_ == 0)
            fieldPtr_->correctImmersedBoundaryConditions();

        correctUnsetBoundaryConditions();

        correct<bcsOfType<patchBoundary>>();
        correct<bcsOfType<parallelBoundary>>();
        correct<bcsOfType<emptyBoundary>>();
    }
}

template<class Type, class MeshType>
template<class Selector, int Degree>
void meshLevel<Type,MeshType>::prepare()
{
    if (fieldPtr_)
    {
        addBoundaryConditions();

        forAll(boundaryConditions_, i)
            if (boundaryConditions_[i].offsetDegree() <= Degree)
                if (Selector::match(boundaryConditions_[i]))
                    boundaryConditions_[i].prepare();
    }
}

template<class Type, class MeshType>
template<class Selector, int Degree>
void meshLevel<Type,MeshType>::evaluate()
{
    if (fieldPtr_)
        forAll(boundaryConditions_, i)
            if (boundaryConditions_[i].offsetDegree() <= Degree)
                if (Selector::match(boundaryConditions_[i]))
                    boundaryConditions_[i].evaluate();
}

template<class Type, class MeshType>
template<class Selector, int Degree>
void meshLevel<Type,MeshType>::correct()
{
    if (fieldPtr_)
    {
        addBoundaryConditions();

        const label nReq = Pstream::nRequests();

        prepare<Selector,Degree>();

        if (Pstream::parRun())
            UPstream::waitRequests(nReq);

        evaluate<Selector,Degree>();
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctUnsetBoundaryConditions()
{
    forAll(*this, d)
    {
        meshDirection<Type,MeshType>& field = listType::operator[](d);

        forAllBlock(this->lvl().boundaries().mask(), i, j, k)
        if (!this->lvl().boundaries().mask()(i,j,k))
        {
            if (labelVector(i,j,k) == unitXYZ)
                continue;

            const labelVector bo(labelVector(i,j,k) - unitXYZ);

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
void meshLevel<Type,MeshType>::correctImmersedBoundaryConditions()
{
    if (fieldPtr_ && fvMsh_.immersedBoundaries<MeshType>().size())
    {
        this->addImmersedBoundaryConditions();

        forAll(fieldPtr_->immersedBoundaryConditions(), i)
        {
            immersedBoundaryCondition<Type,MeshType>& ibc =
                fieldPtr_->immersedBoundaryConditions()[i];

            ibc.evaluate(l_);
        }
    }
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::correctAggData()
{
    const decomposition& decomp = lvl().decomp();

    List<Type> buffer;
    List<labelVector> N(MeshType::numberOfDirections);

    // On agglomerate levels the slaves need to receive from master

    if (decomp.agglomerated())
    {
        const labelBlock& map = decomp.aggProcMap();

        // Prepare buffer. On master we should not account for agglomerate data
        // but on slaves we should.

        label size = 0;
        forAll(N, d)
        {
            N[d] = listType::operator[](d).dataShape(!decomp.aggMaster());
            size += cmptProduct(N[d]);
        }

        buffer.resize(size);

        // Communicate

        if (decomp.aggMaster())
        {
            // Copy data to buffer

            label cursor = 0;
            forAll(*this, d)
                for (int i = 0; i < N[d].x(); i++)
                for (int j = 0; j < N[d].y(); j++)
                for (int k = 0; k < N[d].z(); k++)
                    buffer[cursor++] =
                        listType::operator[](d).B()(i,j,k);

            // Send to slaves

            forAllBlockLinear(map, proc)
                if (proc > 0)
                    UOPstream::write
                    (
                        Pstream::commsTypes::blocking,
                        proc,
                        reinterpret_cast<char*>(buffer.begin()),
                        buffer.byteSize(),
                        0,
                        lvl().comms().agg()
                    );
        }
        else
        {
            // Receive from master

            UIPstream::read
            (
                Pstream::commsTypes::blocking,
                0,
                reinterpret_cast<char*>(buffer.begin()),
                buffer.byteSize(),
                0,
                lvl().comms().agg()
            );

            // Copy data from buffer

            label cursor = 0;
            forAll(*this, d)
                for (int i = 0; i < N[d].x(); i++)
                for (int j = 0; j < N[d].y(); j++)
                for (int k = 0; k < N[d].z(); k++)
                    listType::operator[](d).B()(i,j,k) =
                        buffer[cursor++];
        }
    }

    // On agglomerate parent levels the slaves need to send to master

    if (decomp.aggParent())
    {
        const decomposition& childDecomp = lvl().child().decomp();
        const labelBlock& map = childDecomp.aggProcMap();

        // Prepare buffer

        label size = 0;
        forAll(N, d)
        {
            N[d] = listType::operator[](d).dataShape();
            size += cmptProduct(N[d]);
        }

        buffer.resize(size);

        // Communicate

        if (childDecomp.aggSlave())
        {
            // Copy data to buffer

            label cursor = 0;
            forAll(*this, d)
                for (int i = 0; i < N[d].x(); i++)
                for (int j = 0; j < N[d].y(); j++)
                for (int k = 0; k < N[d].z(); k++)
                    buffer[cursor++] =
                        listType::operator[](d).B()(i,j,k);

            // Send to master

            UOPstream::write
            (
                Pstream::commsTypes::blocking,
                0,
                reinterpret_cast<char*>(buffer.begin()),
                buffer.byteSize(),
                0,
                lvl().child().comms().agg()
            );
        }
        else
        {
            // Receive from slaves

            label proc = 0;
            forAllBlock(map, i, j, k)
            {
                if (proc > 0)
                {
                    const labelVector ijk(i,j,k);

                    UIPstream::read
                    (
                        Pstream::commsTypes::blocking,
                        proc,
                        reinterpret_cast<char*>(buffer.begin()),
                        buffer.byteSize(),
                        0,
                        lvl().child().comms().agg()
                    );

                    // Store

                    label cursor = 0;
                    forAll(*this, d)
                    {
                        meshDirection<Type,MeshType>& dir =
                            listType::operator[](d);

                        // Copy from buffer to block

                        block<Type> B(N[d]);

                        for (int a = 0; a < N[d].x(); a++)
                        for (int b = 0; b < N[d].y(); b++)
                        for (int c = 0; c < N[d].z(); c++)
                            B(a,b,c) =
                                buffer[cursor++];

                        // Start index for the block (offset by ghosts)

                        labelVector Sb(unitXYZ);

                        // Start and end indices for the direction

                        labelVector S(briscola::cmptMultiply(ijk, lvl().N()));
                        labelVector E(S + lvl().N() + MeshType::padding[d]);

                        // Extend into upper ghosts

                        for (int e = 0; e < 3; e++)
                            if (ijk[e] == map.shape()[e]-1)
                                E[e]++;

                        // Extend into lower ghosts

                        for (int e = 0; e < 3; e++)
                        {
                            if (ijk[e] == 0)
                            {
                                S[e]--;
                                Sb[e]--;
                            }
                        }

                        // Copy to direction

                        labelVector abc;
                        for (abc.x() = S.x(); abc.x() < E.x(); abc.x()++)
                        for (abc.y() = S.y(); abc.y() < E.y(); abc.y()++)
                        for (abc.z() = S.z(); abc.z() < E.z(); abc.z()++)
                        {
                            dir(abc) = B(abc-S+Sb);
                        }
                    }
                }

                proc++;
            }
        }
    }
}

template<class Type, class MeshType>
tmp<meshLevel<typename meshLevel<Type,MeshType>::cmptType,MeshType>>
meshLevel<Type,MeshType>::component
(
    const label dir
) const
{
    tmp<meshLevel<cmptType,MeshType>> tL =
        meshLevel<cmptType,MeshType>::New(this->fvMsh_, this->l_);

    forAll(*this, d)
        tL.ref()[d] = listType::operator[](d).component(dir);

    return tL;
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::replace
(
    const label dir,
    const List<cmptType>& values
)
{
    forAll(*this, d)
        listType::operator[](d).replace(dir,values[d]);
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::replace
(
    const label dir,
    const meshLevel<cmptType,MeshType>& L
)
{
    forAll(*this, d)
        listType::operator[](d).replace(dir,L[d]);
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::replace
(
    const label dir,
    const tmp<meshLevel<cmptType,MeshType>>& tL
)
{
    this->replace(dir,tL());

    if (tL.isTmp())
        tL.clear();
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::max(const Type& v)
{
    forAll(*this, d)
        listType::operator[](d).max(v);
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::min(const Type& v)
{
    forAll(*this, d)
        listType::operator[](d).min(v);
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator=(const meshLevel<Type,MeshType>& L)
{
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

    forAll(*this, d)
        listType::operator[](d) = L[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator=
(
    const tmp<meshLevel<Type,MeshType>>& tL
)
{
    #ifdef FULLDEBUG
    checkLevel(tL());
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

    forAll(*this, d)
        listType::operator[](d) += L[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator+=
(
    const tmp<meshLevel<Type,MeshType>>& tL
)
{
    #ifdef FULLDEBUG
    checkLevel(tL());
    #endif

    *this += tL();
    if (tL.isTmp())
        tL.clear();
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator-=(const meshLevel<Type,MeshType>& L)
{
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

    forAll(*this, d)
        listType::operator[](d) -= L[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator-=
(
    const tmp<meshLevel<Type,MeshType>>& tL
)
{
    #ifdef FULLDEBUG
    checkLevel(tL());
    #endif

    *this -= tL();
    if (tL.isTmp())
        tL.clear();
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator*=(const meshLevel<scalar,MeshType>& L)
{
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

    forAll(*this, d)
        listType::operator[](d) *= L[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator*=
(
    const tmp<meshLevel<scalar,MeshType>>& tL
)
{
    #ifdef FULLDEBUG
    checkLevel(tL());
    #endif

    *this *= tL();
    if (tL.isTmp())
        tL.clear();
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator/=(const meshLevel<scalar,MeshType>& L)
{
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

    forAll(*this, d)
        listType::operator[](d) /= L[d];
}

template<class Type, class MeshType>
void meshLevel<Type,MeshType>::operator/=
(
    const tmp<meshLevel<scalar,MeshType>>& tL
)
{
    #ifdef FULLDEBUG
    checkLevel(tL());
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(tL());
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(tL());
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(tL());
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(tL());
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(L);
    #endif

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
    #ifdef FULLDEBUG
    checkLevel(tL());
    #endif

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
