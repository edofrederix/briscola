#include "twoPass.H"
#include "addToRunTimeSelectionTable.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(twoPass, 0);
addToRunTimeSelectionTable(CCL, twoPass, dictionary);

int twoPass::find(label x)
{
    return find(x, unionTable_);
}

int twoPass::find(label x, HashTable<label,label>& table)
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

void twoPass::uniteTags(label x, label y)
{
    label xRoot = find(x);
    label yRoot = find(y);

    if (xRoot != yRoot)
    {
        unionTable_.set
        (
            Foam::max(xRoot, yRoot),
            Foam::min(xRoot, yRoot)
        );
    }
}

void twoPass::parallelReduce()
{
    colocatedLabelField& m = *this;

    m.correctCommsBoundaryConditions();

    List<label> neighbors(0, Zero);

    // Identify equivalent tags in ghost cells

    forAllCells(m,i,j,k)
    {
        if (m(i,j,k))
        {
            neighbors.clear();

            if (m(i-1,j,k))
            {
                neighbors.append(m(i-1,j,k));
            }
            if (m(i,j-1,k))
            {
                neighbors.append(m(i,j-1,k));
            }
            if (m(i,j,k-1))
            {
                neighbors.append(m(i,j,k-1));
            }
            if (m(i+1,j,k))
            {
                neighbors.append(m(i+1,j,k));
            }
            if (m(i,j+1,k))
            {
                neighbors.append(m(i,j+1,k));
            }
            if (m(i,j,k+1))
            {
                neighbors.append(m(i,j,k+1));
            }

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
            IStringStream is(gatheredTables[procI]);
            HashTable<label,label> localTable(is);

            forAll(localTable.toc(), e)
            {
                mergedTable.insert
                (
                    localTable.toc()[e],
                    localTable[localTable.toc()[e]]
                );
            }
        }

        HashTable<label,label> compactTable;
        label nextLabel = 1;

        forAll(mergedTable.toc(), e)
        {
            label root = find(mergedTable.toc()[e], mergedTable);

            if (compactTable.found(root))
            {
                compactTable.insert(mergedTable.toc()[e], compactTable[root]);
            }

            if (!compactTable.found(root))
            {
                compactTable.insert(mergedTable.toc()[e], nextLabel);
                compactTable.insert(root, nextLabel++);
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
        if (m(i,j,k))
        {
            neighbors.clear();

            if (m(i-1,j,k))
            {
                neighbors.append(m(i-1,j,k));
            }
            if (m(i,j-1,k))
            {
                neighbors.append(m(i,j-1,k));
            }
            if (m(i,j,k-1))
            {
                neighbors.append(m(i,j,k-1));
            }

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
    CCL(fvMsh, dict, alpha)
{}

twoPass::twoPass(const twoPass& s)
:
    CCL(s)
{}

twoPass::~twoPass()
{}

}

}

}
