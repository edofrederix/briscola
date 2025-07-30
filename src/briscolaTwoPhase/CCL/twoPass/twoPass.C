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

const label twoPass::procStride = 1e5;

int twoPass::find(label x)
{
    if (!unionTable_.found(x))
    {
        unionTable_.insert(x,x);
    }

    if (unionTable_[x] != x)
    {
        unionTable_.set(x,find(unionTable_[x]));
    }

    return unionTable_[x];
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

twoPass::twoPass
(
    const fvMesh& fvMsh,
    const dictionary& dict,
    const colocatedScalarField& alpha
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

void twoPass::tag()
{
    unionTable_.clear();

    List<label> neighbors(0, Zero);

    label tag = Pstream::myProcNo() * procStride + 1;

    bool firstLoop = true;
    bool changed = true;

    colocatedLabelField& m = *this;

    while (changed)
    {
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

                if (firstLoop)
                {
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
                else
                {
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
        }

        changed = false;

        // Second pass
        forAllCells(m, i,j,k)
        {
            if
            (
                m(i,j,k) &&
                m(i,j,k) != find(m(i,j,k))
            )
            {
                m(i,j,k) = find(m(i,j,k));

                changed = true;
            }
        }

        // Global or reduction, only stop
        // loop if no changes on any processor
        reduce(changed, orOp<bool>());

        // Exchange ghost cells
        m.correctCommsBoundaryConditions();

        if (firstLoop)
            firstLoop = false;
    }
}

}

}

}
