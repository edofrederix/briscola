#include "linearSystemAggregation.H"

#include "SortableList.H"
#include "linearSystem.H"
#include "patchBoundary.H"
#include "parallelBoundary.H"
#include "periodicBoundary.H"
#include "emptyBoundary.H"
#include "PstreamGlobalsLsa.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
linearSystemAggregation<SType,Type,MeshType>::linearSystemAggregation
(
    const linearSystem<SType,Type,MeshType>& sys,
    const label l,
    const label nParts
)
:
    fvMsh_(sys.fvMsh()),
    l_(l),
    nParts_(nParts),
    partStarts_(MeshType::numberOfDirections),
    partEnds_(MeshType::numberOfDirections),
    colNums_(MeshType::numberOfDirections)
{
    if (nParts_ > Pstream::nProcs())
        FatalErrorInFunction
            << "More parts than processors" << endl << abort(FatalError);

    partNum_ = scalar(Pstream::myProcNo()*nParts_)/Pstream::nProcs();
    partMasterNum_ = ceil(scalar(partNum_*Pstream::nProcs())/nParts_);
    nProcsPerPart_ =
        ceil(scalar((partNum_+1)*Pstream::nProcs())/nParts_)
      - partMasterNum_;

    // Pout<< nParts_ << " "
    //     << partNum_ << " "
    //     << partMasterNum_ << " "
    //     << nProcsPerPart_ << " "
    //     << this->master() << " "
    //     << l_ << " "
    //     << d_ << endl;

    // Set LSA communicator

    partCommNum_ =
        PstreamGlobals::lsaSetComms
        (
            nParts_,
            partNum_,
            Pstream::myProcNo() - partMasterNum_
        );

    // Set LSA master communicator

    masterCommNum_ =
        PstreamGlobals::lsaSetComms(nParts_, master() ? 0 : -1, partNum_, true);

    // Set column number lists

    const labelVector* offsets = SType::componentOffsets;

    globalSizes_.setSize(MeshType::numberOfDirections);

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        List<List<FixedList<label,SType::nCsComponents>>>& colNums =
            colNums_[d];

        colNums.setSize(nProcsPerPart_);

        List<FixedList<label,SType::nCsComponents>>& myColNums =
            colNums[Pstream::myProcNo() - partMasterNum_];

        const meshDirection<label,MeshType>& numbers =
            fvMsh_.template metrics<MeshType>().globalCellNumbers()[l_][d];

        myColNums.setSize(numbers.size());

        int c = 0;
        forAllCells(numbers, i, j, k)
        {
            for(int s = 0; s < SType::nCsComponents; s++)
            {
                myColNums[c][s] =
                    numbers(labelVector(i,j,k) + offsets[s]);
            }

            c++;
        }

        // Send/receive column numbers

        if (this->master())
        {
            // Part master, receive from slaves

            for (label proc = 1; proc < nProcsPerPart_; proc++)
            {
                IPstream recv
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    0,
                    0,
                    partCommNum_
                );

                recv >> colNums[proc];
            }
        }
        else
        {
            // Part slave, send to master

            OPstream send
            (
                Pstream::commsTypes::blocking,
                0,
                0,
                0,
                partCommNum_
            );

            send << myColNums;
        }

        // Send/receive global start and end

        if (this->master())
        {
            // Part master, send to slave

            partStarts_[d] = numbers(0,0,0);
            partEnds_[d] = colNums[nProcsPerPart_-1].last()[0] + 1;

            for (label proc = 1; proc < nProcsPerPart_; proc++)
            {
                OPstream send
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    0,
                    0,
                    partCommNum_
                );

                send << partStarts_[d];
                send << partEnds_[d];
            }
        }
        else
        {
            // Part slave, send to master

            IPstream recv
            (
                Pstream::commsTypes::blocking,
                0,
                0,
                0,
                partCommNum_
            );

            recv >> partStarts_[d];
            recv >> partEnds_[d];
        }

        globalSizes_[d] = returnReduce(numbers.size(), sumOp<label>());
    }
}

template<class SType, class Type, class MeshType>
linearSystemAggregation<SType,Type,MeshType>::linearSystemAggregation
(
    const linearSystemAggregation<SType,Type,MeshType>& lsa
)
:
    fvMsh_(lsa.fvMsh_),
    l_(lsa.l_),
    nParts_(lsa.nParts_),
    partNum_(lsa.partNum_),
    partMasterNum_(lsa.partMasterNum_),
    nProcsPerPart_(lsa.nProcsPerPart_),
    globalSizes_(lsa.globalSizes_),
    partStarts_(lsa.partStarts_),
    partEnds_(lsa.partEnds_),
    colNums_(lsa.colNums_),
    partCommNum_(lsa.partCommNum_),
    masterCommNum_(lsa.masterCommNum_)
{}

template<class SType, class Type, class MeshType>
linearSystemAggregation<SType,Type,MeshType>::~linearSystemAggregation()
{}

template<class SType, class Type, class MeshType>
void linearSystemAggregation<SType,Type,MeshType>::rowCoeffs
(
    List<List<SType>>& rows,
    const linearSystem<SType,Type,MeshType>& sys,
    const label d
) const
{
    const label l = this->l_;
    const meshDirection<SType,MeshType>& A = sys.A()[l][d];

    const List<List<FixedList<label,SType::nCsComponents>>>& colNums =
        colNums_[d];

    if (l != l_)
        FatalErrorInFunction
            << "Mesh direction is of invalid level"
            << abort(FatalError);

    const meshDirection<label,MeshType>& numbers =
        fvMsh_.template metrics<MeshType>().globalCellNumbers()[l][d];

    const labelVector* offsets = SType::componentOffsets;

    // Prepare data

    const List<FixedList<label,SType::nCsComponents>>& myColNums =
        colNums_[d][Pstream::myProcNo() - partMasterNum_];

    List<SType> myRows(A.size());

    label c = 0;
    forAllCells(A, i, j, k)
    {
        const labelVector ijk(i,j,k);

        myRows[c] = A(ijk);

        // Remove boundary coefficients

        for (int s = 0; s < SType::nCsComponents; s++)
            if
            (
                numbers(ijk + offsets[s]) < 0
             || numbers(ijk + offsets[s]) >= globalSizes_[d]
            )
                myRows[c][s] = 0;

        // Move redundant indices

        for (int s = 1; s < SType::nCsComponents; s++)
        {
            const int t = findIndex(myColNums[c], myColNums[c][s]);

            if (t < s)
            {
                myRows[c][t] += myRows[c][s];
                myRows[c][s] = 0;
            }
        }

        c++;
    }

    if (this->master())
    {
        // Receive data

        rows.resize(nProcsPerPart_);

        rows[0] = myRows;

        for (label proc = 1; proc < nProcsPerPart_; proc++)
        {
            rows[proc].resize(colNums[proc].size());

            UIPstream::read
            (
                Pstream::commsTypes::blocking,
                proc,
                reinterpret_cast<char*>(rows[proc].begin()),
                rows[proc].byteSize(),
                0,
                partCommNum_
            );
        }
    }
    else
    {
        rows.resize(0);

        // Send data

        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            0,
            reinterpret_cast<char*>(myRows.begin()),
            myRows.byteSize(),
            0,
            partCommNum_
        );
    }
}

template<class SType, class Type, class MeshType>
void linearSystemAggregation<SType,Type,MeshType>::rhsSource
(
    List<Type>& rhs,
    const linearSystem<SType,Type,MeshType>& sys,
    const label d
) const
{
    const label l = this->l_;
    const meshDirection<SType,MeshType>& A = sys.A()[l][d];
    const meshDirection<Type,MeshType>& x = sys.x()[l][d];
    const meshDirection<Type,MeshType>& b = sys.b()[l][d];

    const List<List<FixedList<label,SType::nCsComponents>>>& colNums =
        colNums_[d];

    if (l != l_)
        FatalErrorInFunction
            << "Mesh direction is of invalid level"
            << abort(FatalError);

    const meshDirection<label,MeshType>& numbers =
        fvMsh_.template metrics<MeshType>().globalCellNumbers()[l][d];

    const labelVector* offsets = SType::componentOffsets;

    // Prepare data

    rhs.resize(this->master() ? this->partSize(d) : b.size());

    label c = 0;
    forAllCells(b, i, j, k)
    {
        const labelVector ijk(i,j,k);

        rhs[c] = b(ijk);

        // Add boundary source

        for (int s = 0; s < SType::nCsComponents; s++)
            if
            (
                numbers(ijk + offsets[s]) < 0
             || numbers(ijk + offsets[s]) >= globalSizes_[d]
            )
                rhs[c] -= A(ijk)[s]*x(ijk + offsets[s]);

        c++;
    }

    if (this->master())
    {
        // Receive data

        label c = b.size();

        for (label proc = 1; proc < nProcsPerPart_; proc++)
        {
            const label size = colNums[proc].size();

            UIPstream::read
            (
                Pstream::commsTypes::blocking,
                proc,
                reinterpret_cast<char*>(&rhs[c]),
                size*sizeof(Type),
                0,
                partCommNum_
            );

            c += size;
        }
    }
    else
    {
        // Send

        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            0,
            reinterpret_cast<char*>(rhs.begin()),
            rhs.byteSize(),
            0,
            partCommNum_
        );

        rhs.clear();
    }
}

template<class SType, class Type, class MeshType>
void linearSystemAggregation<SType,Type,MeshType>::compressedRowFormat
(
    scalarList& values,
    labelList& inners,
    labelList& outers,
    const linearSystem<SType,Type,MeshType>& sys,
    const label d,
    const bool sort
) const
{
    const label l = this->l_;

    const label globalStart = partStarts_[d];

    const List<List<FixedList<label,SType::nCsComponents>>>& colNums =
        colNums_[d];

    if (l != l_)
        FatalErrorInFunction
            << "Mesh direction is of invalid level"
            << abort(FatalError);

    List<List<SType>> rows;
    this->rowCoeffs(rows, sys, d);

    if (this->master())
    {
        // Count non-zero values

        int nz = 0;
        forAll(rows, proc)
            forAll(rows[proc], row)
                forAll(rows[proc][row], col)
                    nz += rows[proc][row][col] != 0;

        // Create and set components

        const label size = this->partSize(d);

        values.resize(nz);
        inners.resize(nz);
        outers.resize(size+1);

        outers[0] = 0;

        nz = 0;
        label ri = 0;
        forAll(rows, proc)
        {
            forAll(rows[proc], row)
            {
                forAll(rows[proc][row], col)
                {
                    if (rows[proc][row][col] != 0)
                    {
                        values[nz] = rows[proc][row][col];
                        inners[nz] = colNums[proc][row][col] - globalStart;

                        nz++;
                    }
                }

                outers[++ri] = nz;
            }
        }

        if (sort)
        {
            for (int j = 0; j < size; j++)
            {
                const label s = outers[j];
                const label e = outers[j+1];
                const label n = e - s;

                const SortableList<label> I(SubList<label>(inners, n, s));
                const scalarList V(SubList<scalar>(values, n, s));

                forAll(I, i)
                {
                    inners[s+i] = I[i];
                    values[s+i] = V[I.indices()[i]];
                }
            }
        }
    }
}

template<class SType, class Type, class MeshType>
void linearSystemAggregation<SType,Type,MeshType>::collect
(
    List<Type>& data,
    const meshDirection<Type,MeshType>& f
) const
{
    const label d = f.directionNum();

    const List<List<FixedList<label,SType::nCsComponents>>>& colNums =
        colNums_[d];

    data.resize(this->master() ? this->partSize(d) : f.size());

    label c = 0;
    forAllCells(f, i, j, k)
        data[c++] = f(i,j,k);

    if (this->master())
    {
        // Receive data

        c = f.size();

        for (label proc = 1; proc < nProcsPerPart_; proc++)
        {
            const label size = colNums[proc].size();

            UIPstream::read
            (
                Pstream::commsTypes::blocking,
                proc,
                reinterpret_cast<char*>(&data[c]),
                size*sizeof(Type),
                0,
                partCommNum_
            );

            c += size;
        }
    }
    else
    {
        // Send

        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            0,
            reinterpret_cast<char*>(data.begin()),
            data.byteSize(),
            0,
            partCommNum_
        );

        data.clear();
    }
}

template<class SType, class Type, class MeshType>
void linearSystemAggregation<SType,Type,MeshType>::distribute
(
    meshDirection<Type,MeshType>& f,
    const List<Type>& data
) const
{
    const label d = f.directionNum();

    if (master())
    {
        // Send/copy sol

        label c = 0;
        forAllCells(f, i, j, k)
            f(i,j,k) = data[c++];

        for (label proc = 1; proc < nProcsPerPart_; proc++)
        {
            const label size = this->colNums()[d][proc].size();

            UOPstream::write
            (
                Pstream::commsTypes::blocking,
                proc,
                reinterpret_cast<const char*>(&data[c]),
                size*sizeof(Type),
                0,
                partCommNum_
            );

            c += size;
        }
    }
    else
    {
        // Receive and copy solution

        List<Type> data(f.size());

        UIPstream::read
        (
            Pstream::commsTypes::blocking,
            0,
            reinterpret_cast<char*>(data.begin()),
            data.byteSize(),
            0,
            partCommNum_
        );

        label c = 0;
        forAllCells(f, i, j, k)
            f(i,j,k) = data[c++];
    }
}

// Instantiate

#define INSTANTIATE(STYPE,TYPE,MESHTYPE)                                       \
                                                                               \
    typedef linearSystemAggregation<STYPE,TYPE,MESHTYPE>                       \
        linearSystemAggregation##STYPE##TYPE##MESHTYPE;                        \
                                                                               \
    defineTemplateTypeNameAndDebug                                             \
    (                                                                          \
        linearSystemAggregation##STYPE##TYPE##MESHTYPE,                        \
        0                                                                      \
    )                                                                          \
                                                                               \
    template class linearSystemAggregation<STYPE,TYPE,MESHTYPE>;

INSTANTIATE(diagStencil, scalar, colocated)
INSTANTIATE(diagStencil, scalar, staggered)
INSTANTIATE(diagStencil, vector, colocated)
INSTANTIATE(diagStencil, vector, staggered)

INSTANTIATE(stencil, scalar, colocated)
INSTANTIATE(stencil, scalar, staggered)
INSTANTIATE(stencil, vector, colocated)
INSTANTIATE(stencil, vector, staggered)

#undef INSTANTIATE

}

}

}
