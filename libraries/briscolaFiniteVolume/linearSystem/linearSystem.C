#include "linearSystem.H"
#include "Tuple2.H"

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
        false,
        true
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    singular_(false)
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
        false,
        true
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    singular_(sys.singular_)
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
        false,
        true
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    singular_(tSys->singular_)
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
        false,
        true
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    singular_(sys.singular_)
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
        false,
        true
    ),
    b_
    (
        IOobject::groupName("b", x_.name()),
        fvMsh_,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    singular_(tSys->singular_)
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

    tSys.clear();
}

template<class SType, class Type, class MeshType>
linearSystem<SType,Type,MeshType>::~linearSystem()
{}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::residual
(
    meshField<Type, MeshType>& res
) const
{
    forAll(res, l)
        this->residual(res[l]);
}

template<class SType, class Type, class MeshType>
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

    this->residual(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::residual
(
    meshLevel<Type,MeshType>& res
) const
{
    forAll(res, d)
        this->residual(res[d]);
}

template<class SType, class Type, class MeshType>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::residual(const label l) const
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
) const
{
    const label l = res.levelNum();
    const label d = res.directionNum();

    Amul(res, this->A()[l][d], x_[l][d]);

    res *= -1.0;
    res += this->b()[l][d];
}

template<class SType, class Type, class MeshType>
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

    this->residual(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::evaluate
(
    meshField<Type, MeshType>& res
) const
{
    forAll(res, l)
        this->evaluate(res[l]);
}

template<class SType, class Type, class MeshType>
tmp<meshField<Type, MeshType>>
linearSystem<SType,Type,MeshType>::evaluate() const
{
    tmp<meshField<Type,MeshType>> tRes
    (
        new meshField<Type,MeshType>
        (
            "evaluate",
            fvMsh_
        )
    );

    this->evaluate(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::evaluate
(
    meshLevel<Type,MeshType>& res
) const
{
    forAll(res, d)
        this->evaluate(res[d]);
}

template<class SType, class Type, class MeshType>
tmp<meshLevel<Type,MeshType>>
linearSystem<SType,Type,MeshType>::evaluate(const label l) const
{
    tmp<meshLevel<Type,MeshType>> tRes
    (
        new meshLevel<Type,MeshType>(fvMsh_,l)
    );

    this->evaluate(tRes.ref());

    return tRes;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::evaluate
(
    meshDirection<Type,MeshType>& res
) const
{
    const label l = res.levelNum();
    const label d = res.directionNum();

    const meshDirection<scalar,MeshType>& cv =
        fvMsh_.template metrics<MeshType>().cellVolumes()[l][d];

    Amul(res, this->A()[l][d], x_[l][d]);

    res -= this->b()[l][d];
    res /= cv;
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
    forAll(A_, l)
        this->eliminateGhosts(l);
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::eliminateGhosts(const label l)
{
    forAll(x_.boundaryConditions(), i)
        x_.boundaryConditions()[i].eliminateGhosts(*this, l);
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator=
(
    const linearSystem<SType,Type,MeshType>& sys
)
{
    this->A() = sys.A();
    this->b() = sys.b();
    singular_ = sys.singular();
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
    singular_ = sys.singular() ? singular_ : false;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    *this += tSys();
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
    singular_ = sys.singular() ? singular_ : false;
}

template<class SType, class Type, class MeshType>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const tmp<linearSystem<SType,Type,MeshType>>& tSys
)
{
    *this -= tSys();
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
    singular_ = sys.singular();
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator=
(
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys
)
{
    *this = tSys();
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
    singular_ = sys.singular() ? singular_ : false;
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator+=
(
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys
)
{
    *this += tSys();
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
    singular_ = sys.singular() ? singular_ : false;
}

template<class SType, class Type, class MeshType>
template<class SType2>
void linearSystem<SType,Type,MeshType>::operator-=
(
    const tmp<linearSystem<SType2,Type,MeshType>>& tSys
)
{
    *this -= tSys();
    tSys.clear();
}

template<class SType, class Type, class MeshType>
void writeToFile
(
    linearSystem<SType,Type,MeshType>& sys,
    OFstream& file,
    const label l
)
{
    bool makeShallowAgain = false;

    if (l > 0 && sys.b().shallow())
    {
        sys.b().makeDeep();

        autoPtr<restrictionScheme<Type,MeshType>> reScheme
        (
            restrictionScheme<Type,MeshType>::New("linear", sys.fvMsh())
        );

        for (int m = 1; m <= l; m++)
        for (int d = 0; d < MeshType::numberOfDirections; d++)
        {
            reScheme->restrict(sys.b()[m][d], sys.b()[m-1][d], true);
        }

        makeShallowAgain = true;
    }

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        OFstream * file2Ptr;

        if (d == 0)
        {
            file2Ptr = &file;
        }
        else
        {
            file2Ptr = new OFstream(file.name() + "_" + Foam::name(d));
        }

        OFstream& file2 = *file2Ptr;

        const meshDirection<SType,MeshType>& A = sys.A()[l][d];
        const meshDirection<Type,MeshType>& b = sys.b()[l][d];

        const labelVector N = A.N();

        List<labelVector> shapes(Pstream::nProcs());
        shapes[Pstream::myProcNo()] = N;

        Pstream::gatherList(shapes);

        List<List<Tuple2<stencil,Type>>> data(Pstream::nProcs());

        data[Pstream::myProcNo()].setSize(cmptProduct(N));

        forAllCells(A, i, j, k)
        {
            const label l = i*N.y()*N.z() + j*N.z() + k;

            data[Pstream::myProcNo()][l].first() = A(i,j,k);
            data[Pstream::myProcNo()][l].second() = b(i,j,k);
        }

        Pstream::gatherList(data);

        if (Pstream::master())
        {
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
                    const stencil& Ai = data[proc][i].first();
                    const Type& bi = data[proc][i].second();

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
                            file2 << Ai[0] << " ";
                        }
                        else
                        {
                            bool found = false;
                            for (int f = 0; f < 6; f++)
                            {
                                if (neigh[f] == j)
                                {
                                    file2 << Ai[f+1] << " ";
                                    found = true;
                                    break;
                                }
                            }

                            if (!found)
                                file2 << 0.0 << " ";
                        }
                    }

                    file2 << bi << endl;
                }
            }
        }

        if (d > 0)
            free(file2Ptr);
    }

    if (makeShallowAgain)
        sys.b().makeShallow();
}

}

}

}
