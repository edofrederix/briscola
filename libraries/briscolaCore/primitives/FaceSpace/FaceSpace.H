#ifndef FaceSpace_H
#define FaceSpace_H

#include "CellSpace.H"
#include "briscolaUtils.H"
#include "pTransform.H"
#include "stencil.H"
#include "LowerFaceSpace.H"

namespace Foam
{

template<class Cmpt>
class FaceSpace
:
    public CellSpace<FaceSpace<Cmpt>, Cmpt, 6>
{

public:

    // Component type
    typedef Cmpt cmpt;

    // Component of component type
    typedef scalar cmptCmpt;

    // Face space type
    typedef FaceSpace<Cmpt> fsType;

    // Cell space type
    typedef CellSpace<FaceSpace<Cmpt>, Cmpt, 6> csType;

    // Component labeling enumeration
    enum components { L, R, B, T, A, F };

    // Construct null
    inline FaceSpace();

    // Construct initialised to zero
    inline FaceSpace(const Foam::zero);

    // Construct given cell space of the same size
    template<class Cmpt2>
    inline FaceSpace(const CellSpace<FaceSpace<Cmpt2>, Cmpt2, 6>&);

    // Construct from six components
    inline FaceSpace
    (
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&
    );

    // Construct from lower and upper vectors
    inline FaceSpace(const Vector<Cmpt>&, const Vector<Cmpt>&);

    // Construct from Istream
    inline FaceSpace(Istream&);

    // Access to left, right, bottom, top, aft and fore element

    inline const Cmpt& left() const { return this->v_[L]; }
    inline const Cmpt& right() const { return this->v_[R]; }
    inline const Cmpt& bottom() const { return this->v_[B]; }
    inline const Cmpt& top() const { return this->v_[T]; }
    inline const Cmpt& aft() const { return this->v_[A]; }
    inline const Cmpt& fore() const { return this->v_[F]; }

    inline Cmpt& left()  { return this->v_[L]; }
    inline Cmpt& right()  { return this->v_[R]; }
    inline Cmpt& bottom()  { return this->v_[B]; }
    inline Cmpt& top()  { return this->v_[T]; }
    inline Cmpt& aft()  { return this->v_[A]; }
    inline Cmpt& fore()  { return this->v_[F]; }

    // Return lower elements as vector
    inline Vector<Cmpt> lower() const
    {
        return Vector<Cmpt>(left(), bottom(), aft());
    }

    // Return upper elements as vector
    inline Vector<Cmpt> upper() const
    {
        return Vector<Cmpt>(right(), top(), fore());
    }

    // Conversion to lower face space
    inline operator LowerFaceSpace<Cmpt>() const;

    // Conversion to stencil (only implemented for Cmpt = scalar)
    inline operator stencil() const;

    // Explicitly use operators from cell space
    using csType::operator=;
    using csType::operator+=;
    using csType::operator-=;
    using csType::operator*=;
    using csType::operator/=;
};

// FaceSpace partial specialization of CellSpaceCast

template<class Cmpt1, class Cmpt2, class Cmpt3>
class CellSpaceCast<FaceSpace<Cmpt1>, FaceSpace<Cmpt2>, Cmpt3>
{
public:
    typedef FaceSpace<Cmpt3> type;
};

template<class Cmpt1, class Cmpt2, class Cmpt3>
class CellSpaceCast<LowerFaceSpace<Cmpt1>, FaceSpace<Cmpt2>, Cmpt3>
{
public:
    typedef LowerFaceSpace<Cmpt3> type;
};

template<class Cmpt1, class Cmpt2, class Cmpt3>
class CellSpaceCast<FaceSpace<Cmpt1>, LowerFaceSpace<Cmpt2>, Cmpt3>
{
public:
    typedef LowerFaceSpace<Cmpt3> type;
};

// FaceSpace partial specialization of pTransform

using briscola::faceNumber;
using briscola::faceOffsets;

template<class Cmpt>
class pTransform<FaceSpace<Cmpt>>
{
    FaceSpace<Cmpt> p_;

public:

    inline pTransform(const labelTensor& T, const FaceSpace<Cmpt>& h)
    {
        if (T == briscola::eye)
        {
            p_ = h;
        }
        else
        {
            for (int i = 0; i < 6; i++)
                p_[max(faceNumber(T & faceOffsets[i]), 0)] = h[i];
        }
    }

    inline operator FaceSpace<Cmpt>()
    {
        return p_;
    }
};

}

#include "FaceSpaceProducts.H"
#include "FaceSpaceI.H"

#endif
