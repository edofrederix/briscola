#ifndef partPoints_H
#define partPoints_H

// Wrapper around vectorBlock, which adds ghost cell padding useful for periodic
// or parallel boundaries. The base vector block is kept private.

#include "labelVector.H"
#include "block.H"

namespace Foam
{

namespace briscola
{

class partPoints
:
    private vectorBlock
{
public:

    // Null constructor
    partPoints();

    // Destructor
    ~partPoints();

    // Clear the point level
    void clear();

    // Set the size of the level point block from the number of cells, i.e.,
    // excluding ghost cells
    void setSizeFromCells(const labelVector&);

    // Access to the points, padded by the ghost point layer so that (0,0,0)
    // corresponds to the left-bottom-aft point of the first cell

    inline vector& operator()(const int i, const int j, const int k)
    {
        return vectorBlock::operator()(i+1,j+1,k+1);
    }

    inline const vector& operator()(const int i, const int j, const int k) const
    {
        return vectorBlock::operator()(i+1,j+1,k+1);
    }

    inline vector& operator()(const labelVector& ijk)
    {
        return vectorBlock::operator()(ijk + unitXYZ);
    }

    inline const vector& operator()(const labelVector& ijk) const
    {
        return vectorBlock::operator()(ijk + unitXYZ);
    }

    // Scalar interpolation, useful for non-colocated meshes

    inline vector interp(const scalar i, const scalar j, const scalar k) const
    {
        return vectorBlock::interp(i+1,j+1,k+1);
    }

    inline vector interp(const vector ijk) const
    {
        return vectorBlock::interp(ijk+vector(unitXYZ));
    }

    // Linear element access operator will not work efficiently due to
    // ghost cells, so let's prevent its use

    inline vector& operator()(const int i)
    {
        NotImplemented;
        return vectorBlock::operator()(i);
    }

    inline const vector& operator()(const int i) const
    {
        NotImplemented;
        return vectorBlock::operator()(i);
    }

    // Return the first dimension without ghost points
    inline int l() const { return vectorBlock::l()-2; }

    // Return the second dimension
    inline int m() const { return vectorBlock::m()-2; }

    // Return the third dimension
    inline int n() const { return vectorBlock::n()-2; }

    // Return the number of elements in the block without ghost points
    inline int size() const { return l()*m()*n(); }

    // Return the size of the block without ghost points
    inline labelVector shape() const
    {
        return labelVector(l(),m(),n());
    }

    // Alias for shape()
    inline labelVector N() const { return shape(); }

    // Return the points of given cell
    inline vertexVector cellPoints
    (
        const label i,
        const label j,
        const label k
    ) const
    {
        return vertexVector
        (
            this->operator()(i  ,j  ,k  ),
            this->operator()(i+1,j  ,k  ),
            this->operator()(i  ,j+1,k  ),
            this->operator()(i+1,j+1,k  ),
            this->operator()(i  ,j  ,k+1),
            this->operator()(i+1,j  ,k+1),
            this->operator()(i  ,j+1,k+1),
            this->operator()(i+1,j+1,k+1)
        );
    }

    inline vertexVector cellPoints(const labelVector& ijk) const
    {
        return cellPoints(ijk.x(), ijk.y(), ijk.z());
    }

    // Compute the local coordinates of a point in the given cell, calculating
    // the weights of the tri-linear interpolation function in the cell.
    inline vector cellCoordinates
    (
        const vector& p,
        const label i,
        const label j,
        const label k,
        const bool insideOnly = false
    ) const
    {
        return interpolationWeights(p,this->cellPoints(i,j,k),insideOnly);
    }

    inline vector cellCoordinates
    (
        const vector& p,
        const labelVector& ijk,
        const bool insideOnly = false
    ) const
    {
        return interpolationWeights(p,this->cellPoints(ijk),insideOnly);
    }

    // Check if a point is inside the given cell, which is true if the local
    // coordinates are on the unit cube [0,1]^3
    inline bool pointInCell
    (
        const vector& p,
        const label i,
        const label j,
        const label k
    ) const
    {
        return this->cellCoordinates(p,i,j,k,true) != -vector::one;
    }

    inline bool pointInCell(const vector& p, const labelVector& ijk) const
    {
        return this->cellCoordinates(p,ijk,true) != -vector::one;
    }
};

}

}

#endif
