#ifndef EdgeSpace_H
#define EdgeSpace_H

#include "CellSpace.H"
#include "blockUtils.H"
#include "pTransform.H"
#include "stencil.H"

namespace Foam
{

template<class Cmpt>
class EdgeSpace
:
    public CellSpace<EdgeSpace<Cmpt>, Cmpt, 12>
{

public:

    // Component type
    typedef Cmpt cmpt;

    // Component of component type
    typedef scalar cmptCmpt;

    // Edge space type
    typedef EdgeSpace<Cmpt> fsType;

    // Cell space type
    typedef CellSpace<EdgeSpace<Cmpt>, Cmpt, 12> csType;

    // Component labeling enumeration
    enum components { X0, X1, X2, X3, Y0, Y1, Y2, Y3, Z0, Z1, Z2, Z3 };

    // Construct null
    inline EdgeSpace();

    // Construct initialised to zero
    inline EdgeSpace(const Foam::zero);

    // Construct uniform
    inline EdgeSpace(const Cmpt&);

    // Construct given cell space of the same size
    template<class Cmpt2>
    inline EdgeSpace(const CellSpace<EdgeSpace<Cmpt2>, Cmpt2, 12>&);

    // Construct from six components
    inline EdgeSpace
    (
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&,
        const Cmpt&
    );

    // Construct from Istream
    inline EdgeSpace(Istream&);

    // Access to edge elements

    inline const Cmpt& x0() const { return this->v_[X0]; }
    inline const Cmpt& x1() const { return this->v_[X1]; }
    inline const Cmpt& x2() const { return this->v_[X2]; }
    inline const Cmpt& x3() const { return this->v_[X3]; }

    inline const Cmpt& y0() const { return this->v_[Y0]; }
    inline const Cmpt& y1() const { return this->v_[Y1]; }
    inline const Cmpt& y2() const { return this->v_[Y2]; }
    inline const Cmpt& y3() const { return this->v_[Y3]; }

    inline const Cmpt& z0() const { return this->v_[Z0]; }
    inline const Cmpt& z1() const { return this->v_[Z1]; }
    inline const Cmpt& z2() const { return this->v_[Z2]; }
    inline const Cmpt& z3() const { return this->v_[Z3]; }

    inline Cmpt& x0() { return this->v_[X0]; }
    inline Cmpt& x1() { return this->v_[X1]; }
    inline Cmpt& x2() { return this->v_[X2]; }
    inline Cmpt& x3() { return this->v_[X3]; }

    inline Cmpt& y0() { return this->v_[Y0]; }
    inline Cmpt& y1() { return this->v_[Y1]; }
    inline Cmpt& y2() { return this->v_[Y2]; }
    inline Cmpt& y3() { return this->v_[Y3]; }

    inline Cmpt& z0() { return this->v_[Z0]; }
    inline Cmpt& z1() { return this->v_[Z1]; }
    inline Cmpt& z2() { return this->v_[Z2]; }
    inline Cmpt& z3() { return this->v_[Z3]; }
};

// EdgeSpace partial specialization of CellSpaceCast

template<class Cmpt1, class Cmpt2, class Cmpt3>
class CellSpaceCast<EdgeSpace<Cmpt1>, EdgeSpace<Cmpt2>, Cmpt3>
{
public:
    typedef EdgeSpace<Cmpt3> type;
};

// EdgeSpace partial specialization of pTransform

template<class Cmpt>
class pTransform<EdgeSpace<Cmpt>>
{
    EdgeSpace<Cmpt> p_;

public:

    inline pTransform(const labelTensor& T, const EdgeSpace<Cmpt>& h)
    {
        if (T == briscola::eye)
        {
            p_ = h;
        }
        else
        {
            const labelVector x0(T & briscola::edgeOffset0);
            const labelVector x1(T & briscola::edgeOffset1);

            const labelVector y0(T & briscola::edgeOffset4);
            const labelVector y1(T & briscola::edgeOffset5);

            const labelVector z0(T & briscola::edgeOffset8);
            const labelVector z1(T & briscola::edgeOffset9);

            p_[briscola::edgeNumber( x0)] = h.x0();
            p_[briscola::edgeNumber( x1)] = h.x1();
            p_[briscola::edgeNumber(-x1)] = h.x2();
            p_[briscola::edgeNumber(-x0)] = h.x3();

            p_[briscola::edgeNumber( y0)] = h.y0();
            p_[briscola::edgeNumber( y1)] = h.y1();
            p_[briscola::edgeNumber(-y1)] = h.y2();
            p_[briscola::edgeNumber(-y0)] = h.y3();

            p_[briscola::edgeNumber( z0)] = h.z0();
            p_[briscola::edgeNumber( z1)] = h.z1();
            p_[briscola::edgeNumber(-z1)] = h.z2();
            p_[briscola::edgeNumber(-z0)] = h.z3();
        }
    }

    inline operator EdgeSpace<Cmpt>()
    {
        return p_;
    }
};

}

#include "EdgeSpaceProducts.H"
#include "EdgeSpaceI.H"

#endif