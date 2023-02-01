#ifndef block_H
#define block_H

#include "zero.H"
#include "autoPtr.H"
#include "List.H"
#include "tmp.H"
#include "refCount.H"
#include "labelVector.H"
#include "labelTensor.H"
#include "products.H"

#include "blockFwd.H"
#include "blockUtils.H"

namespace Foam
{

namespace briscola
{

template<class Type> class block;

template<class Type>
Istream& operator>>(Istream&, block<Type>&);

template<class Type>
Ostream& operator<<(Ostream&, const block<Type>&);

template<class Type>
class block
:
    public refCount
{
    // Dimensions of the block
    int l_;
    int m_;
    int n_;

    // Row pointers
    Type* __restrict__ v_;

    // Transformation of the block w.r.t. its construction
    labelTensor T_;

public:

    // Component type
    typedef typename pTraits<Type>::cmptType cmptType;

    // Null constructor
    inline block()
    :
        l_(0),
        m_(0),
        n_(0),
        v_(nullptr),
        T_(eye)
    {}

    // Construct with given dimensions
    block(const int, const int, const int);

    // Construct with given dimensions and initialize as zero
    block(const int, const int, const int, const zero);

    // Construct with given dimensions and initialize with value
    block(const int, const int, const int, const Type&);

    // Construct with given shape
    block(const labelVector&);

    // Construct with given shape and initialize as zero
    block(const labelVector&, const zero);

    // Construct with given shape and initialize with value
    block(const labelVector&, const Type&);

    // Construct with given dimensions as copy of a list
    block(const int, const int, const int, const List<Type>&);

    // Construct with given dimensions as copy of an array (doesn't check
    // bounds)
    block(const int, const int, const int, const Type*);

    // Construct with given shape as copy of a list
    block(const labelVector&, const List<Type>&);

    // Construct with given shape as copy of an array (doesn't check bounds)
    block(const labelVector&, const Type*);

    // Copy constructors (preserves rotation)
    block(const block<Type>&);
    block(const block<Type>&, const zero&);
    block(const block<Type>&, const Type&);

    // Copy constructor from tmp and delete argument (preserves rotation)
    block(const tmp<block<Type>>&);
    block(const tmp<block<Type>>&, const zero&);
    block(const tmp<block<Type>>&, const Type&);

    // Construct from Istream
    block(Istream& is);

    // Clone
    inline autoPtr<block<Type>> clone() const
    {
        return autoPtr<block<Type>>(new block<Type>(*this));
    }

    // Destructor
    ~block();

    // Allocate the storage for the element vector
    inline void allocate();

    // Clear the block, i.e., set sizes to zero
    inline void clear();

    // Re-allocate the storage for the element vector
    inline void reAllocate(int, int, int);

    // Re-allocate the storage for the element vector
    inline void reAllocate(const labelVector&);

    // Transfer the data of the argument block into this block and annul the
    // argument block. Rotation is left as-is
    inline void transferData(block<Type>&);

    // Return the first dimension
    inline int l() const { return l_; }

    // Return the second dimension
    inline int m() const { return m_; }

    // Return the third dimension
    inline int n() const { return n_; }

    // Size of the block in bytes
    inline std::streamsize byteSize() const
    {
        return size()*sizeof(Type);
    }

    // Start of the array
    inline Type* begin() { return v_; }

    // Const start of the array
    inline const Type* cbegin() const { return v_; }

    // End of the array
    inline Type* end() { return &v_[size()]; }

    // Const end of the array
    inline const Type* cend() const { return &v_[size()]; }

    // Return acces to the transformation tensor
    inline labelTensor& T() { return T_; }
    inline const labelTensor& T() const { return T_; }

    #ifdef FULLDEBUG

    // Check if given linear i index is within bounds
    inline void checkLineari(const int i) const
    {
        if (i < 0 || i >= l_*m_*n_)
        {
            FatalErrorInFunction
                << "Linear index i = " << i << " out of range 0 ... " << l_*m_*n_-1
                << abort(FatalError);
        }
    }

    // Check if given i index is within bounds
    inline void checki(const int i) const
    {
        if (i < 0 || i >= l_)
        {
            FatalErrorInFunction
                << "Index i = " << i << " out of range 0 ... " << l_-1
                << abort(FatalError);
        }
    }

    // Check if given j index is within bounds
    inline void checkj(const int j) const
    {
        if (j < 0 || j >= m_)
        {
            FatalErrorInFunction
                << "Index j = " << j << " out of range 0 ... " << m_-1
                << abort(FatalError);
        }
    }

    // Check if given k index is within bounds
    inline void checkk(const int k) const
    {
        if (k < 0 || k >= n_)
        {
            FatalErrorInFunction
                << "Index k = " << k << " out of range 0 ... " << n_-1
                << abort(FatalError);
        }
    }

    #endif

    // Return the number of elements in the block
    inline int size() const { return l_*m_*n_; }

    // Return the shape of the block
    inline labelVector shape() const
    {
        return labelVector(l_,m_,n_);
    }

    // Return the shape of the block
    inline labelVector N() const
    {
        return shape();
    }

    // Return the shape of the block as a diagonal tensor
    inline labelTensor shapeTensor() const
    {
        return labelTensor(l_,0,0,0,m_,0,0,0,n_);
    }

    // Const access to the element vector
    inline const Type* v() const { return v_; }

    // Access to the element vector
    inline Type* v() { return v_; }

    // Resize the block preserving the elements
    void setSize(const labelVector);

    // Resize the block preserving the elements
    void setSize(const int, const int, const int);

    // Prepend a surface in given direction with given thickness and uniform
    // value
    void prepend(const int, const int, const Type&);

    // Append a surface in given direction with given thickness and uniform
    // value
    void append(const int, const int, const Type&);

    // Transform the block with given matrix transformation (rotation,
    // reflection or permutation)
    void transform(const labelTensor);

    // Reflect array in given axis and return transformation tensor
    labelTensor reflect(const label);

    // Permute axes and return transformation tensor
    labelTensor permute(const label, const label);

    // Rotate the block in counter-clockwise direction about given axis and by
    // given number of steps and return the transformation tensor
    labelTensor rotate(const label, const label);

    // Squeeze the block (move non-unity dimensions forward)
    void squeeze();

    // Slice block perpendicular to the x-axis at given index. Optionally
    // squeeze the block. Negative index allowed
    tmp<block<Type>> slicex(const label, const bool = false) const;

    // Slice block perpendicular to the y-axis at given index. Optionally
    // squeeze the block. Negative index allowed
    tmp<block<Type>> slicey(const label, const bool = false) const;

    // Slice block perpendicular to the z-axis at given index. Optionally
    // squeeze the block. Negative index allowed
    tmp<block<Type>> slicez(const label, const bool = false) const;

    // Slice block at given index and perpendicular to given direction.
    // Optionally squeeze the block. Negative index allowed
    tmp<block<Type>> slice
    (
        const label,
        const label,
        const bool = false
    ) const;

    // Slice block at given edge number. Optionally squeeze the block
    tmp<block<Type>> sliceEdge(const label, const bool = false) const;

    // Slice block at given vertex number
    tmp<block<Type>> sliceVertex(const label) const;

    // Find the position of given value in the block. Linear search, so not a
    // good idea to use on large blocks
    labelVector find(const Type&) const;

    // Compute linear index from components
    inline int lin
    (
        const label& i,
        const label& j,
        const label& k
    ) const
    {
        return i*m_*n_ + j*n_ + k;
    }

    // Linear index element access operator
    inline Type& operator()(const int i)
    {
        #ifdef FULLDEBUG

        checkLineari(i);

        #endif

        return v_[i];
    }

    // Linear index const element access operator
    inline const Type& operator()(const int i) const
    {
        #ifdef FULLDEBUG

        checkLineari(i);

        #endif

        return v_[i];
    }

    // (i,j,k) element access operator
    inline Type& operator()
    (
        const int i,
        const int j,
        const int k
    )
    {
        #ifdef FULLDEBUG

        checki(i);
        checkj(j);
        checkk(k);

        #endif

        return v_[lin(i,j,k)];
    }

    // (i,j,k) const element access operator
    inline const Type& operator()
    (
        const int i,
        const int j,
        const int k
    ) const
    {
        #ifdef FULLDEBUG

        checki(i);
        checkj(j);
        checkk(k);

        #endif

        return v_[lin(i,j,k)];
    }

    // (i,j,k) element access operator
    inline Type& operator()(const labelVector& ijk)
    {
        #ifdef FULLDEBUG

        checki(ijk.x());
        checkj(ijk.y());
        checkk(ijk.z());

        #endif

        return v_[lin(ijk.x(),ijk.y(),ijk.z())];
    }

    // (i,j,k) const element access operator
    inline const Type& operator()(const labelVector& ijk) const
    {
        #ifdef FULLDEBUG

        checki(ijk.x());
        checkj(ijk.y());
        checkk(ijk.z());

        #endif

        return v_[lin(ijk.x(),ijk.y(),ijk.z())];
    }

    // Interpolate data to a point shifted from given indices. The magnitude of
    // each component of the shift vector must be in between or equal to -1 and
    // 1

    Type interp
    (
        const int i,
        const int j,
        const int k,
        const vector& shift
    ) const;

    Type interp
    (
        const labelVector& ijk,
        const vector& shift
    ) const;

    // Member operators

    inline void operator=(const block<Type>&);
    inline void operator=(const tmp<block<Type>>&);
    inline void operator=(const Type&);
    inline void operator=(const zero);

    inline void operator+=(const block<Type>&);
    inline void operator+=(const tmp<block<Type>>&);

    inline void operator-=(const block<Type>&);
    inline void operator-=(const tmp<block<Type>>&);

    inline void operator*=(const block<scalar>&);
    inline void operator*=(const tmp<block<scalar>>&);

    inline void operator/=(const block<scalar>&);
    inline void operator/=(const tmp<block<scalar>>&);

    inline void operator+=(const Type&);
    inline void operator-=(const Type&);

    inline void operator*=(const scalar&);
    inline void operator/=(const scalar&);

    // Templated functions are needed for types which can be cast to another
    // type

    template<class Type2> inline void operator=(const block<Type2>&);
    template<class Type2> inline void operator=(const tmp<block<Type2>>&);
    template<class Type2> inline void operator+=(const block<Type2>&);
    template<class Type2> inline void operator+=(const tmp<block<Type2>>&);
    template<class Type2> inline void operator-=(const block<Type2>&);
    template<class Type2> inline void operator-=(const tmp<block<Type2>>&);
    template<class Type2> inline void operator=(const Type2&);
    template<class Type2> inline void operator+=(const Type2&);
    template<class Type2> inline void operator-=(const Type2&);

    friend Istream& operator>><Type>(Istream&, block<Type>&);
    friend Ostream& operator<<<Type>(Ostream&, const block<Type>&);
};

}

}

#include "blockFunctions.H"
#include "blockI.H"

#ifdef NoRepository
    #include "block.C"
#endif

#endif
