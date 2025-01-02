#include "linearSystemAggregation.H"

#include "SortableList.H"
#include "linearSystem.H"
#include "domainBoundary.H"
#include "parallelBoundary.H"
#include "periodicBoundary.H"
#include "emptyBoundary.H"

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
    singular_(const_cast<linearSystem<SType,Type,MeshType>&>(sys).singular()),
    globalStarts_(MeshType::numberOfDirections),
    globalEnds_(MeshType::numberOfDirections),
    colNums_(MeshType::numberOfDirections)
{
    if (nParts_ > Pstream::nProcs())
        FatalErrorInFunction
            << "More parts than processors" << endl << abort(FatalError);

    myPartNum_ = scalar(Pstream::myProcNo()*nParts_)/Pstream::nProcs();
    myPartMasterNum_ = ceil(scalar(myPartNum_*Pstream::nProcs())/nParts_);
    nProcsPerPart_ =
        ceil(scalar((myPartNum_+1)*Pstream::nProcs())/nParts_)
      - myPartMasterNum_;

    // Pout<< nParts_ << " "
    //     << myPartNum_ << " "
    //     << myPartMasterNum_ << " "
    //     << nProcsPerPart_ << " "
    //     << this->master() << " "
    //     << l_ << " "
    //     << d_ << endl;

    for (int d = 0; d < MeshType::numberOfDirections; d++)
    {
        List<List<FixedList<label,FullSType::nComponents>>>& colNums =
            colNums_[d];

        colNums.setSize(nProcsPerPart_);

        List<FixedList<label,FullSType::nComponents>>& myColNums =
            colNums[Pstream::myProcNo() - myPartMasterNum_];

        const meshDirection<label,MeshType>& numbers =
            fvMsh_.template metrics<MeshType>().globalCellNumbers()[l_][d];

        myColNums.setSize(numbers.size());

        int c = -1;
        forAllCells(numbers, i, j, k)
            for(int s = 0; s < FullSType::nComponents; s++)
                myColNums[s == 0 ? ++c : c][s] =
                    numbers
                    (
                        labelVector(i,j,k)
                      + FullSType::componentOffsets[s]
                    );

        // Send/receive column numbers

        if (this->master())
        {
            // Part master, receive from slaves

            for
            (
                label proc = myPartMasterNum_ + 1;
                proc < myPartMasterNum_ + nProcsPerPart_;
                proc++
            )
            {
                IPstream recv(Pstream::commsTypes::blocking, proc);
                recv >> colNums[proc - myPartMasterNum_];
            }
        }
        else
        {
            // Part slave, send to master

            OPstream send(Pstream::commsTypes::blocking, myPartMasterNum_);
            send << myColNums;
        }

        // Send/receive global start and end

        if (this->master())
        {
            // Part master, send to slave

            globalStarts_[d] = numbers(0,0,0);
            globalEnds_[d] = colNums[nProcsPerPart_-1].last()[0] + 1;

            for
            (
                label proc = myPartMasterNum_ + 1;
                proc < myPartMasterNum_ + nProcsPerPart_;
                proc++
            )
            {
                OPstream send(Pstream::commsTypes::blocking, proc);
                send << globalStarts_[d];
                send << globalEnds_[d];
            }
        }
        else
        {
            // Part slave, send to master

            IPstream recv(Pstream::commsTypes::blocking, myPartMasterNum_);
            recv >> globalStarts_[d];
            recv >> globalEnds_[d];
        }
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
    singular_(lsa.singular_),
    myPartNum_(lsa.myPartNum_),
    myPartMasterNum_(lsa.myPartMasterNum_),
    nProcsPerPart_(lsa.nProcsPerPart_),
    globalStarts_(lsa.globalStarts_),
    globalEnds_(lsa.globalEnds_),
    colNums_(lsa.colNums_)
{}

template<class SType, class Type, class MeshType>
linearSystemAggregation<SType,Type,MeshType>::~linearSystemAggregation()
{}

template<class SType, class Type, class MeshType>
void linearSystemAggregation<SType,Type,MeshType>::rowCoeffs
(
    List<List<FullSType>>& rows,
    const meshDirection<SType,MeshType>& A
) const
{
    const label l = A.levelNum();
    const label d = A.directionNum();

    const label globalStart = globalStarts_[d];
    const label globalEnd = globalEnds_[d];

    const List<List<FixedList<label,FullSType::nComponents>>>& colNums =
        colNums_[d];

    if (l != l_)
        FatalErrorInFunction
            << "Mesh direction is of invalid level"
            << abort(FatalError);

    const meshDirection<label,MeshType>& numbers =
        fvMsh_.template metrics<MeshType>().globalCellNumbers()[l][d];

    const labelVector* offsets = FullSType::componentOffsets;

    // Prepare data

    const List<FixedList<label,FullSType::nComponents>>& myColNums =
        colNums_[d][Pstream::myProcNo() - myPartMasterNum_];

    List<FullSType> myRows(A.size());

    label c = 0;
    forAllCells(A, i, j, k)
    {
        const labelVector ijk(i,j,k);

        myRows[c] = fullStencil(A,ijk);

        // Remove boundary coefficients

        for (int s = 0; s < FullSType::nComponents; s++)
            if
            (
                numbers(ijk + offsets[s]) < globalStart
             || numbers(ijk + offsets[s]) >= globalEnd
            )
                myRows[c][s] = 0;

        // Add constraint if needed

        if (singular_[d] && myPartMasterNum_ == 0)
            for (int s = 0; s < FullSType::nComponents; s++)
                if (numbers(ijk + FullSType::componentOffsets[s]) == 0)
                    myRows[c][s] = 0;

        // Move redundant indices

        for (int s = 1; s < FullSType::nComponents; s++)
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

    // Add constraint if needed

    if (singular_[d] && Pstream::myProcNo() == 0)
        myRows[0][0] = 1;

    if (this->master())
    {
        // Receive data

        rows.resize(nProcsPerPart_);

        for
        (
            label proc = myPartMasterNum_;
            proc < myPartMasterNum_ + nProcsPerPart_;
            proc++
        )
        {
            const label i = proc-myPartMasterNum_;
            rows[i].resize(colNums[i].size());

            if (i == 0)
            {
                rows[i] = myRows;
            }
            else
            {
                UIPstream::read
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    reinterpret_cast<char*>(rows[i].begin()),
                    rows[i].byteSize(),
                    0,
                    UPstream::worldComm
                );
            }
        }
    }
    else
    {
        rows.resize(0);

        // Send data

        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            myPartMasterNum_,
            reinterpret_cast<char*>(myRows.begin()),
            myRows.byteSize(),
            0,
            UPstream::worldComm
        );
    }
}

template<class SType, class Type, class MeshType>
void linearSystemAggregation<SType,Type,MeshType>::rhsSource
(
    List<Type>& rhs,
    const meshDirection<SType,MeshType>& A,
    const meshDirection<Type,MeshType>& x,
    const meshDirection<Type,MeshType>& b
) const
{
    const label l = A.levelNum();
    const label d = A.directionNum();

    const label globalStart = globalStarts_[d];
    const label globalEnd = globalEnds_[d];

    const List<List<FixedList<label,FullSType::nComponents>>>& colNums =
        colNums_[d];

    if (l != l_)
        FatalErrorInFunction
            << "Mesh direction is of invalid level"
            << abort(FatalError);

    const meshDirection<label,MeshType>& numbers =
        fvMsh_.template metrics<MeshType>().globalCellNumbers()[l][d];

    const labelVector* offsets = FullSType::componentOffsets;

    // Prepare data

    List<Type> myRhs(b.size());

    label c = 0;
    forAllCells(b, i, j, k)
    {
        const labelVector ijk(i,j,k);

        myRhs[c] = b(ijk);

        // Add boundary source

        for (int s = 0; s < FullSType::nComponents; s++)
            if
            (
                numbers(ijk + offsets[s]) < globalStart
             || numbers(ijk + offsets[s]) >= globalEnd
            )
                myRhs[c] -=
                    fullStencil(A,ijk)[s]
                  * x(ijk + FullSType::componentOffsets[s]);

        c++;
    }

    // Add constraint if needed

    if (singular_[d] && Pstream::myProcNo() == 0)
        myRhs[0] = Zero;

    if (this->master())
    {
        // Receive data

        rhs.resize(this->size(d));

        label c = 0;

        for
        (
            label proc = myPartMasterNum_;
            proc < myPartMasterNum_ + nProcsPerPart_;
            proc++
        )
        {
            const label i = proc - myPartMasterNum_;
            const label size = colNums[i].size();

            if (i == 0)
            {
                forAll(myRhs, j)
                    rhs[j] = myRhs[j];
            }
            else
            {
                UIPstream::read
                (
                    Pstream::commsTypes::blocking,
                    proc,
                    reinterpret_cast<char*>(&rhs[c]),
                    size*sizeof(Type),
                    0,
                    UPstream::worldComm
                );
            }

            c += size;
        }
    }
    else
    {
        rhs.resize(0);

        // Send

        UOPstream::write
        (
            Pstream::commsTypes::blocking,
            myPartMasterNum_,
            reinterpret_cast<char*>(myRhs.begin()),
            myRhs.byteSize(),
            0,
            UPstream::worldComm
        );
    }
}

template<class SType, class Type, class MeshType>
void linearSystemAggregation<SType,Type,MeshType>::compressedRowFormat
(
    scalarList& values,
    labelList& inners,
    labelList& outers,
    const meshDirection<SType,MeshType>& A,
    const bool sort
) const
{
    const label l = A.levelNum();
    const label d = A.directionNum();

    const label globalStart = globalStarts_[d];

    const List<List<FixedList<label,FullSType::nComponents>>>& colNums =
        colNums_[d];

    if (l != l_)
        FatalErrorInFunction
            << "Mesh direction is of invalid level"
            << abort(FatalError);

    List<List<FullSType>> rows;
    this->rowCoeffs(rows, A);

    if (this->master())
    {
        // Count non-zero values

        int nz = 0;
        forAll(rows, proc)
            forAll(rows[proc], row)
                forAll(rows[proc][row], col)
                    nz += rows[proc][row][col] != 0;

        // Create and set components

        const label size = this->size(d);

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
void linearSystemAggregation<SType,Type,MeshType>::correctBoundaryConditions
(
    meshLevel<Type,MeshType>& x
) const
{
    meshField<Type,MeshType>& f = x.mshField();

    const label l = x.levelNum();

    if (l != l_)
        FatalErrorInFunction
            << "Mesh level is of invalid level"
            << abort(FatalError);

    f.addBoundaryConditions();

    // First non-eliminated domain boundaries

    forAll(f.boundaryConditions(), i)
    {
        boundaryCondition<Type,MeshType>& bc = f.boundaryConditions()[i];
        const boundary& b = bc.mshBoundary();

        if (!bc.eliminated() && b.castable<domainBoundary>())
            bc.evaluate(l);
    }

    // Next the parallel/periodic boundaries between the aggregated matrices

    const label nReq = Pstream::nRequests();

    forAll(f.boundaryConditions(), i)
    {
        boundaryCondition<Type,MeshType>& bc = f.boundaryConditions()[i];
        const boundary& b = bc.mshBoundary();

        if (b.castable<parallelBoundary>())
        {
            const label procNum = b.cast<parallelBoundary>().neighborProcNum();

            if
            (
                procNum < myPartMasterNum_
             || procNum >= myPartMasterNum_ + nProcsPerPart_
             || b.castable<periodicBoundary>()
            )
            {
                bc.prepare(l);
            }
        }
    }

    if (Pstream::parRun())
    {
        Pstream::waitRequests(nReq);
    }

    forAll(f.boundaryConditions(), i)
    {
        boundaryCondition<Type,MeshType>& bc = f.boundaryConditions()[i];
        const boundary& b = bc.mshBoundary();

        if (b.castable<parallelBoundary>())
        {
            const label procNum = b.cast<parallelBoundary>().neighborProcNum();

            if
            (
                procNum < myPartMasterNum_
             || procNum >= myPartMasterNum_ + nProcsPerPart_
             || b.castable<periodicBoundary>()
            )
            {
                bc.evaluate(l);
            }
        }
    }

    // Finally the other non-eliminated boundaries

    forAll(f.boundaryConditions(), i)
    {
        boundaryCondition<Type,MeshType>& bc = f.boundaryConditions()[i];
        const boundary& b = bc.mshBoundary();

        if
        (
            !bc.eliminated()
         && !b.castable<domainBoundary>()
         && !b.castable<parallelBoundary>()
        )
        {
            bc.evaluate(l);
        }
    }
}

}

}

}
