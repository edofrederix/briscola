#ifndef VertexSpace_H
#define VertexSpace_H

#include "CellSpace.H"
#include "utils.H"
#include "pTransform.H"
#include "stencil.H"

namespace Foam
{

template<class Cmpt>
class VertexSpace
:
    public CellSpace<VertexSpace<Cmpt>, Cmpt, 8>
{

public:

    // Component type
    typedef Cmpt cmpt;

    // Component of component type
    typedef scalar cmptCmpt;

    // Edge space type
    typedef VertexSpace<Cmpt> fsType;

    // Cell space type
    typedef CellSpace<VertexSpace<Cmpt>, Cmpt, 8> csType;

    // Component labeling enumeration
    enum components { LBA, RBA, LTA, RTA, LBF, RBF, LTF, RTF };

    // Construct null
    inline VertexSpace();

    // Construct initialised to zero
    inline VertexSpace(const Foam::zero);

    // Construct uniform
    inline VertexSpace(const Cmpt&);

    // Construct given cell space of the same size
    template<class Cmpt2>
    inline VertexSpace(const CellSpace<VertexSpace<Cmpt2>, Cmpt2, 8>&);

    // Construct from six components
    inline VertexSpace
    (
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
    inline VertexSpace(Istream&);

    // Access to vertex elements

    inline const Cmpt& lba() const { return this->v_[LBA]; }
    inline const Cmpt& rba() const { return this->v_[RBA]; }
    inline const Cmpt& lta() const { return this->v_[LTA]; }
    inline const Cmpt& rta() const { return this->v_[RTA]; }
    inline const Cmpt& lbf() const { return this->v_[LBF]; }
    inline const Cmpt& rbf() const { return this->v_[RBF]; }
    inline const Cmpt& ltf() const { return this->v_[LTF]; }
    inline const Cmpt& rtf() const { return this->v_[RTF]; }

    inline Cmpt& lba() { return this->v_[LBA]; }
    inline Cmpt& rba() { return this->v_[RBA]; }
    inline Cmpt& lta() { return this->v_[LTA]; }
    inline Cmpt& rta() { return this->v_[RTA]; }
    inline Cmpt& lbf() { return this->v_[LBF]; }
    inline Cmpt& rbf() { return this->v_[RBF]; }
    inline Cmpt& ltf() { return this->v_[LTF]; }
    inline Cmpt& rtf() { return this->v_[RTF]; }
};

// VertexSpace partial specialization of CellSpaceCast

template<class Cmpt1, class Cmpt2, class Cmpt3>
class CellSpaceCast<VertexSpace<Cmpt1>, VertexSpace<Cmpt2>, Cmpt3>
{
public:
    typedef VertexSpace<Cmpt3> type;
};

// VertexSpace partial specialization of pTransform

using briscola::vertexNumber;
using briscola::vertexOffsets;

template<class Cmpt>
class pTransform<VertexSpace<Cmpt>>
{
    VertexSpace<Cmpt> p_;

public:

    inline pTransform(const labelTensor& T, const VertexSpace<Cmpt>& h)
    {
        if (T == briscola::eye)
        {
            p_ = h;
        }
        else
        {
            for (int i = 0; i < 8; i++)
                p_[max(vertexNumber(T & vertexOffsets[i]), 0)] = h[i];
        }
    }

    inline operator VertexSpace<Cmpt>()
    {
        return p_;
    }
};

}

#include "VertexSpaceProducts.H"
#include "VertexSpaceI.H"

#endif