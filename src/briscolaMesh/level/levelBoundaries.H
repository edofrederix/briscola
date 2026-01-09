#ifndef levelBoundaries_H
#define levelBoundaries_H

#include "FastPtrList.H"
#include "boundary.H"
#include "labelBlock.H"
#include "brickLink.H"

namespace Foam
{

namespace briscola
{

class level;

class levelBoundaries
:
    public FastPtrList<boundary>
{
    // Level to which these boundaries belong
    const level& lvl_;

    // Mask to keep track of set boundaries
    labelBlock mask_;

    // Mask to keep track of parallel/periodic boundaries
    labelBlock parallelMask_;

    // Add a boundary to the list of boundaries
    void add
    (
        const word name,
        const word type,
        const labelTensor T,
        const labelVector offset,
        const vector mappedOffset = vector::zero,
        const labelVector neighborOffset = zeroXYZ,
        const label neighborProcNum = -1
    );

    // Generate brick-internal boundaries
    void generateBrickInternalBoundaries();

    // Generate brick-external boundary
    void generateBrickExternalBoundaries();

    // Generate link boundaries from brick link
    void generateLinkBoundaries(const brickLink&);

    // Generate patch boundaries
    void generatePatchBoundaries();

    // Set boundary mask
    void setMask();

    // Reorder boundaries in the order:
    //
    //  1) patch/mapped boundaries,
    //  2) parallel/periodic faces,
    //  3) parallel/periodic edges,
    //  4) parallel/periodic vertices and
    //  5) empty boundaries.
    //
    //  This order is important, so that parallel boundary conditions can copy
    //  patch/mapped boundaries and empty boundaries can copy both patch/mapped
    //  and parallel/periodic boundaries.
    void reorder();

    // Set communication tags in parallel and periodic boundaries. Tags avoid
    // inconsistent communication.
    void setCommTags();

public:

    // Construct from level
    levelBoundaries(const level&);

    // Copy constructor
    levelBoundaries(const levelBoundaries&);

    // Copy constructor with alternative level reference
    levelBoundaries(const levelBoundaries&, const level&);

    // Destructor
    ~levelBoundaries();

    // Level
    inline const level& lvl() const { return lvl_; }

    // Return the boundary mask
    inline const labelBlock& mask() const
    {
        return mask_;
    }

    // Return the parallel boundary mask
    inline const labelBlock& parallelMask() const
    {
        return parallelMask_;
    }

    // Find boundary with given offset vector
    inline boundary& find(const labelVector& bo)
    {
        forAll(*this, i)
            if (this->operator[](i).offset() == bo)
                return this->operator[](i);

        FatalErrorInFunction
            << "Invalid boundary offset" << endl << abort(FatalError);

        return this->operator[](0);
    }

    // Const reference find with given offset vector
    inline const boundary& find(const labelVector& bo) const
    {
        return const_cast<levelBoundaries*>(this)->find(bo);
    }

    // Find boundary with given name
    inline boundary& find(const word name)
    {
        forAll(*this, i)
            if (this->operator[](i).name() == name)
                return this->operator[](i);

        FatalErrorInFunction
            << "Invalid boundary name" << endl << abort(FatalError);

        return this->operator[](0);
    }

    // Const reference find with given name
    inline const boundary& find(const word name) const
    {
        return const_cast<levelBoundaries*>(this)->find(name);
    }
};

}

}

#endif
