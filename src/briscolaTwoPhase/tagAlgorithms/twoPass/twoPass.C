#include "twoPass.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(twoPass, 0);
addToRunTimeSelectionTable(tagAlgorithm, twoPass, dictionary);

label twoPass::find(const label x)
{
    return find(x, unionTable_);
}

label twoPass::find(const label x, HashTable<label,label>& table)
{
    if (!table.found(x))
    {
        table.insert(x,x);
    }

    if (table[x] != x)
    {
        table.set(x, find(table[x], table));
    }

    return table[x];
}

void twoPass::uniteTags(const label x, const label y)
{
    uniteTags(unionTable_,x,y);
}

void twoPass::uniteTags
(
    HashTable<label,label>& table,
    const label x,
    const label y
)
{
    label xRoot = find(x,table);
    label yRoot = find(y,table);

    if (xRoot != yRoot)
    {
        table.set
        (
            Foam::max(xRoot, yRoot),
            Foam::min(xRoot, yRoot)
        );
    }
}

void twoPass::parallelReduce()
{
    colocatedLabelField& m = *this;

    m.correct<bcsOfType<parallelBoundary>>();

    List<label> neighbors(0, Zero);

    // Identify equivalent tags in ghost cells

    forAllCells(m,i,j,k)
    {
        const labelVector ijk(i,j,k);

        if (m(ijk))
        {
            // Initialize to maximum number of neighbors
            neighbors.setSize(6, Zero);

            label cursor = 0;

            for (int f = 0; f < 6; f++)
                if (m(ijk + faceOffsets[f]))
                    neighbors[cursor++] = m(ijk + faceOffsets[f]);

            // Resize to real number of neighbors
            neighbors.resize(cursor);

            if (neighbors.size() != 0)
            {
                label minLabel = Foam::min(neighbors);

                forAll(neighbors, n)
                {
                    uniteTags(minLabel,neighbors[n]);
                }
            }
        }
    }

    // Serialize HashTable's

    OStringStream os;

    os << unionTable_;

    string localSerialized = os.str();

    // Gather tables

    List<string> gatheredTables(Pstream::nProcs());

    gatheredTables[Pstream::myProcNo()] = localSerialized;

    Pstream::gatherList(gatheredTables);

    // Merge tables on master

    HashTable<label,label> mergedTable;

    string mergedStr;

    if (Pstream::master())
    {
        forAll(gatheredTables, procI)
        {
            if (gatheredTables[procI].empty()) continue;

            IStringStream is(gatheredTables[procI]);
            HashTable<label,label> localTable(is);

            forAllConstIter(labelTable, localTable, iter)
            {
                uniteTags
                (
                    mergedTable,
                    iter.key(),
                    *iter
                );
            }
        }

        HashTable<label,label> compactTable(mergedTable.size());
        label nextLabel = 1;

        forAllConstIter(labelTable, mergedTable, iter)
        {
            const label key = iter.key();

            label root = find(*iter, mergedTable);

            if (!compactTable.found(root))
            {
                compactTable.insert(root, nextLabel);
                compactTable.insert(key, nextLabel);
                nextLabel++;
            }
            else
            {
                compactTable.insert(key, compactTable[root]);
            }
        }

        this->n_ = nextLabel - 1;

        OStringStream mergedStream;

        mergedStream << compactTable;

        mergedStr = mergedStream.str();
    }

    // Scatter tables to all processors

    Pstream::scatter(mergedStr);
    Pstream::scatter(this->n_);

    IStringStream mergedStrStream(mergedStr);

    HashTable<label,label> globalMergedTable(mergedStrStream);

    unionTable_.transfer(globalMergedTable);

    // Resolve unified tags again

    forAllCells(m, i,j,k)
    {
        if (m(i,j,k))
        {
            m(i,j,k) = unionTable_[m(i,j,k)];
        }
    }
}

void twoPass::tag()
{
    unionTable_.clear();

    labelList neighbors(0, Zero);

    label tag =
        fvMsh().metrics<colocated>().globalCellNumbers()(0,0,0) + 1;

    colocatedLabelField& m = *this;

    // First pass
    forAllCells(m, i,j,k)
    {
        const labelVector ijk(i,j,k);

        if (m(ijk))
        {
            // Initialize to maximum number of lower neighbors
            neighbors.setSize(3, Zero);

            label cursor = 0;

            for (int f = 0; f < 3; f++)
                if (m(ijk + lowerFaceOffsets[f]))
                    neighbors[cursor++] = m(ijk + lowerFaceOffsets[f]);

            // Resize to real number of neighbors
            neighbors.resize(cursor);

            if (neighbors.size() == 0)
            {
                m(i,j,k) = tag++;
            }
            else
            {
                label minLabel = Foam::min(neighbors);

                m(i,j,k) = minLabel;

                forAll(neighbors, n)
                {
                    uniteTags(minLabel,neighbors[n]);
                }
            }
        }
    }

    // Second pass
    forAllCells(m, i,j,k)
    {
        if (m(i,j,k))
        {
            m(i,j,k) = find(m(i,j,k));
        }
    }

    // Parallel reduction
    parallelReduce();
}

twoPass::twoPass
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    colocatedScalarField& alpha
)
:
    tagAlgorithm(fvMsh, dict, alpha)
{}

twoPass::twoPass(const twoPass& s)
:
    tagAlgorithm(s)
{}

twoPass::~twoPass()
{}

}

}

}
