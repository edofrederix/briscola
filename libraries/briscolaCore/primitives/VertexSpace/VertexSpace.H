#ifndef VertexSpace_H
#define VertexSpace_H

#include "CellSpace.H"
#include "blockUtils.H"
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
    enum components { V0, V1, V2, V3, V4, V5, V6, V7 };

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

    inline const Cmpt& v0() const { return this->v_[V0]; }
    inline const Cmpt& v1() const { return this->v_[V1]; }
    inline const Cmpt& v2() const { return this->v_[V2]; }
    inline const Cmpt& v3() const { return this->v_[V3]; }
    inline const Cmpt& v4() const { return this->v_[V4]; }
    inline const Cmpt& v5() const { return this->v_[V5]; }
    inline const Cmpt& v6() const { return this->v_[V6]; }
    inline const Cmpt& v7() const { return this->v_[V7]; }

    inline Cmpt& v0() { return this->v_[V0]; }
    inline Cmpt& v1() { return this->v_[V1]; }
    inline Cmpt& v2() { return this->v_[V2]; }
    inline Cmpt& v3() { return this->v_[V3]; }
    inline Cmpt& v4() { return this->v_[V4]; }
    inline Cmpt& v5() { return this->v_[V5]; }
    inline Cmpt& v6() { return this->v_[V6]; }
    inline Cmpt& v7() { return this->v_[V7]; }
};

// VertexSpace partial specialization of CellSpaceCast

template<class Cmpt1, class Cmpt2, class Cmpt3>
class CellSpaceCast<VertexSpace<Cmpt1>, VertexSpace<Cmpt2>, Cmpt3>
{
public:
    typedef VertexSpace<Cmpt3> type;
};

// VertexSpace partial specialization of pTransform

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
            const labelVector v0(T & briscola::vertexOffset0);
            const labelVector v1(T & briscola::vertexOffset1);
            const labelVector v2(T & briscola::vertexOffset2);
            const labelVector v3(T & briscola::vertexOffset3);
            const labelVector v4(T & briscola::vertexOffset4);
            const labelVector v5(T & briscola::vertexOffset5);
            const labelVector v6(T & briscola::vertexOffset6);
            const labelVector v7(T & briscola::vertexOffset7);

            p_[briscola::vertexNumber(v0)] = h.v0();
            p_[briscola::vertexNumber(v1)] = h.v1();
            p_[briscola::vertexNumber(v2)] = h.v2();
            p_[briscola::vertexNumber(v3)] = h.v3();
            p_[briscola::vertexNumber(v4)] = h.v4();
            p_[briscola::vertexNumber(v5)] = h.v5();
            p_[briscola::vertexNumber(v6)] = h.v6();
            p_[briscola::vertexNumber(v7)] = h.v7();
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