#ifndef LowerFaceSpace_H
#define LowerFaceSpace_H

#include "CellSpace.H"
#include "briscolaUtils.H"
#include "pTransform.H"
#include "stencil.H"

namespace Foam
{

template<class Cmpt>
class LowerFaceSpace
:
    public CellSpace<LowerFaceSpace<Cmpt>, Cmpt, 3>
{

public:

    // Component type
    typedef Cmpt cmpt;

    // Component of component type
    typedef scalar cmptCmpt;

    // Face space type
    typedef LowerFaceSpace<Cmpt> fsType;

    // Cell space type
    typedef CellSpace<LowerFaceSpace<Cmpt>, Cmpt, 3> csType;

    // Component labeling enumeration
    enum components { L, B, A };

    // Construct null
    inline LowerFaceSpace();

    // Construct initialised to zero
    inline LowerFaceSpace(const Foam::zero);

    // Construct uniform
    inline LowerFaceSpace(const Cmpt&);

    // Construct given cell space of the same size
    template<class Cmpt2>
    inline LowerFaceSpace(const CellSpace<LowerFaceSpace<Cmpt2>, Cmpt2, 3>&);

    // Construct from three components
    inline LowerFaceSpace
    (
        const Cmpt&,
        const Cmpt&,
        const Cmpt&
    );

    // Construct from lower vector
    inline LowerFaceSpace(const Vector<Cmpt>&);

    // Construct from Istream
    inline LowerFaceSpace(Istream&);

    // Access to left, bottom and aft element

    inline const Cmpt& left() const { return this->v_[L]; }
    inline const Cmpt& bottom() const { return this->v_[B]; }
    inline const Cmpt& aft() const { return this->v_[A]; }

    inline Cmpt& left()  { return this->v_[L]; }
    inline Cmpt& bottom()  { return this->v_[B]; }
    inline Cmpt& aft()  { return this->v_[A]; }

    // Return lower elements as vector
    inline Vector<Cmpt> lower() const
    {
        return Vector<Cmpt>(left(), bottom(), aft());
    }
};

}

#include "LowerFaceSpaceProducts.H"
#include "LowerFaceSpaceI.H"

#endif
