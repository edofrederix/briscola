#ifndef tessellation_H
#define tessellation_H

#include "tetra.H"
#include "hexa.H"
#include "piped.H"

namespace Foam
{

namespace briscola
{

/**

@brief Tessellation geometric object

A tessellation is defined as a list of tetrahedrals (tets). For tets, geometric
operations such as volume computation and truncation are defined. Decomposing
other geometrical objects into a tessellation thus allows to extend such
operations. For example, the volume of a general hexahedral (hex) can be
computed from the sum of the volumes of the tets that decompose this hex.

*/

class tessellation
:
    public List<tetra>
{
public:

    // Null construct
    inline tessellation()
    :
        List<tetra>()
    {}

    // Copy constructor
    inline tessellation(const tessellation& T)
    :
        List<tetra>(T)
    {}

    // Construct from a single tet
    inline tessellation(const tetra& t0)
    :
        List<tetra>(1)
    {
        this->operator[](0) = t0;
    }

    // Construct from a hex
    inline tessellation(const hexa& h)
    :
        List<tetra>(numberOfTets)
    {
        for (label i = 0; i < numberOfTets; i++)
            this->operator[](i) =
                tetra
                (
                    h.v()[tetDecomp[i][0]],
                    h.v()[tetDecomp[i][1]],
                    h.v()[tetDecomp[i][2]],
                    h.v()[tetDecomp[i][3]]
                );
    }

    // Construct from a piped
    inline tessellation(const piped& p)
    :
        List<tetra>(tessellation(hexa(p)))
    {}

    // Construct from two tets
    inline tessellation(const tetra& t0, const tetra& t1)
    :
        List<tetra>(2)
    {
        this->operator[](0) = t0;
        this->operator[](1) = t1;
    }

    // Construct from three tets
    inline tessellation(const tetra& t0, const tetra& t1, const tetra& t2)
    :
        List<tetra>(3)
    {
        this->operator[](0) = t0;
        this->operator[](1) = t1;
        this->operator[](2) = t2;
    }

    // Construct from a list of tets
    inline tessellation(const List<tetra>& tets)
    :
        List<tetra>(tets)
    {}

    // Assignment
    inline void operator=(const tessellation& t)
    {
        static_cast<List<tetra>&>(*this) = t;
    }

    // Compute the center as the volume average of the tet centers
    inline vector center() const
    {
        vector c = Zero;
        scalar v = Zero;

        forAll(*this, i)
        {
            const scalar w = this->operator[](i).volume();

            c += w*this->operator[](i).center();
            v += w;
        }

        return v > 1e-12 ? c/v : Zero;
    }

    // Compute the volume
    inline scalar volume() const
    {
        scalar V = 0;

        forAll(*this, i)
            V += this->operator[](i).volume();

        return V;
    }

    // Check if tesselation is valid
    inline bool valid() const
    {
        forAll(*this, i)
            if (!this->operator[](i).valid())
                return false;

        return true;
    }

    // Compute a new truncated tessellation from this tessellation and a plane
    // defined by dot(x,n) + C = 0
    inline tessellation truncate(const vector& n, const scalar C) const
    {
        tessellation T;
        forAll(*this, i)
            T.append(this->operator[](i).truncate(n,C));

        return T;
    }

    // General intersection function
    template<class Type>
    inline tessellation intersect(const Type& t) const
    {
        return intersect(tessellation(t));
    }

    // Compute intersection between this tessellation and a tet
    inline tessellation intersect(const tetra& t) const
    {
        tessellation R;
        forAll(*this, i)
            R.append(this->operator[](i).intersect(t));

        return R;
    }

    // Compute intersection between this tessellation and another
    inline tessellation intersect(const tessellation& T) const
    {
        tessellation R;
        forAll(*this, i)
            forAll(T, j)
                R.append(this->operator[](i).intersect(T[j]));

        return R;
    }

    // Volume of this tessellation truncated by a place
    inline scalar truncationVolume(const vector& n, const scalar C) const
    {
        return truncate(n,C).volume();
    }

    // Volume of an intersection
    template<class Type>
    inline scalar intersectionVolume(const Type& t) const
    {
        return intersect(t).volume();
    }
};

}

}

#endif
