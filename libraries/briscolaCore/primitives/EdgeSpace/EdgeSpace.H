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
    enum components { BA, TA, BF, TF, LA, RA, LF, RF, LB, RB, LT, RT };

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

    inline const Cmpt& ba() const { return this->v_[BA]; }
    inline const Cmpt& ta() const { return this->v_[TA]; }
    inline const Cmpt& bf() const { return this->v_[BF]; }
    inline const Cmpt& tf() const { return this->v_[TF]; }
    inline const Cmpt& la() const { return this->v_[LA]; }
    inline const Cmpt& ra() const { return this->v_[RA]; }
    inline const Cmpt& lf() const { return this->v_[LF]; }
    inline const Cmpt& rf() const { return this->v_[RF]; }
    inline const Cmpt& lb() const { return this->v_[LB]; }
    inline const Cmpt& rb() const { return this->v_[RB]; }
    inline const Cmpt& lt() const { return this->v_[LT]; }
    inline const Cmpt& rt() const { return this->v_[RT]; }

    inline Cmpt& ba() { return this->v_[BA]; }
    inline Cmpt& ta() { return this->v_[TA]; }
    inline Cmpt& bf() { return this->v_[BF]; }
    inline Cmpt& tf() { return this->v_[TF]; }
    inline Cmpt& la() { return this->v_[LA]; }
    inline Cmpt& ra() { return this->v_[RA]; }
    inline Cmpt& lf() { return this->v_[LF]; }
    inline Cmpt& rf() { return this->v_[RF]; }
    inline Cmpt& lb() { return this->v_[LB]; }
    inline Cmpt& rb() { return this->v_[RB]; }
    inline Cmpt& lt() { return this->v_[LT]; }
    inline Cmpt& rt() { return this->v_[RT]; }
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
            const labelVector ba(T & briscola::edgeOffset0);
            const labelVector ta(T & briscola::edgeOffset1);
            const labelVector bf(T & briscola::edgeOffset2);
            const labelVector tf(T & briscola::edgeOffset3);
            const labelVector la(T & briscola::edgeOffset4);
            const labelVector ra(T & briscola::edgeOffset5);
            const labelVector lf(T & briscola::edgeOffset6);
            const labelVector rf(T & briscola::edgeOffset7);
            const labelVector lb(T & briscola::edgeOffset8);
            const labelVector rb(T & briscola::edgeOffset9);
            const labelVector lt(T & briscola::edgeOffset10);
            const labelVector rt(T & briscola::edgeOffset11);

            p_[briscola::edgeNumber(ba)] = h.ba();
            p_[briscola::edgeNumber(ta)] = h.ta();
            p_[briscola::edgeNumber(bf)] = h.bf();
            p_[briscola::edgeNumber(tf)] = h.tf();
            p_[briscola::edgeNumber(la)] = h.la();
            p_[briscola::edgeNumber(ra)] = h.ra();
            p_[briscola::edgeNumber(lf)] = h.lf();
            p_[briscola::edgeNumber(rf)] = h.rf();
            p_[briscola::edgeNumber(lb)] = h.lb();
            p_[briscola::edgeNumber(rb)] = h.rb();
            p_[briscola::edgeNumber(lt)] = h.lt();
            p_[briscola::edgeNumber(rt)] = h.rt();
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