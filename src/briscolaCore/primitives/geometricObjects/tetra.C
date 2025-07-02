#include "tetra.H"

#include "hexa.H"
#include "piped.H"
#include "tessellation.H"
#include "labelList.H"

namespace Foam
{

namespace briscola
{

tessellation tetra::truncate(const vector& n, const scalar C) const
{
    // No need to avoid truncation of a tet with zero volume

    if (volume() < 1e-12)
        return tessellation();

    scalarList d(4);
    labelList inside(4), outside(4);

    label j = 0, k = 0;
    forAll(d, i)
    {
        d[i] = (v_[i] & n) + C;

        if (d[i] >= 0.0)
        {
            inside[j++] = i;
        }
        else
        {
            outside[k++] = i;
        }
    }

    inside.resize(j);
    outside.resize(k);

    tessellation T;

    // Return empty tessellation if no points are inside. Otherwise, before
    // adding tets, check if they have non-zero volume.

    if (inside.size() == 1)
    {
        const vector p0 = v_[inside[0]];

        const vector p00 = interp(d, inside[0], outside[0]);
        const vector p01 = interp(d, inside[0], outside[1]);
        const vector p02 = interp(d, inside[0], outside[2]);

        const tetra t =
            inside[0]%2 == 0
          ? tetra(p02,p01,p00,p0)
          : tetra(p00,p01,p02,p0);

        if (t.volume() >= 1e-12)
            T.append(t);
    }
    else if (inside.size() == 2)
    {
        const vector p0 = v_[inside[0]];
        const vector p1 = v_[inside[1]];

        const vector p00 = interp(d, inside[0], outside[0]);
        const vector p01 = interp(d, inside[0], outside[1]);
        const vector p10 = interp(d, inside[1], outside[0]);
        const vector p11 = interp(d, inside[1], outside[1]);

        List<tetra> tets(3);

        if ((inside[0] + inside[1])%2 == 0)
        {
            tets[0] = tetra(p11, p00, p01, p0);
            tets[1] = tetra(p00, p11, p10, p1);
            tets[2] = tetra(p00, p11, p1,  p0);
        }
        else
        {
            tets[0] = tetra(p00, p11, p01, p0);
            tets[1] = tetra(p11, p00, p10, p1);
            tets[2] = tetra(p00, p11, p0,  p1);
        }

        forAll(tets, i)
            if (tets[i].volume() >= 1e-12)
                T.append(tets[i]);
    }
    else if (inside.size() == 3)
    {
        const vector p0 = v_[inside[0]];
        const vector p1 = v_[inside[1]];
        const vector p2 = v_[inside[2]];

        const vector p00 = interp(d, inside[0], outside[0]);
        const vector p10 = interp(d, inside[1], outside[0]);
        const vector p20 = interp(d, inside[2], outside[0]);

        List<tetra> tets(3);

        if (outside[0]%2 == 0)
        {
            tets[0] = tetra(p2, p1,  p0,  p00);
            tets[1] = tetra(p2, p1,  p00, p10);
            tets[2] = tetra(p2, p10, p00, p20);

        }
        else
        {
            tets[0] = tetra(p0,  p1,  p2, p00);
            tets[1] = tetra(p00, p1,  p2, p10);
            tets[2] = tetra(p00, p10, p2, p20);
        }

        forAll(tets, i)
            if (tets[i].volume() >= 1e-12)
                T.append(tets[i]);
    }
    else if (inside.size() == 4)
    {
        T.append(tetra(*this));
    }

    return T;
}

tessellation tetra::intersect(const tetra& t) const
{
    tessellation T(*this);

    // Compute the truncation of this tet with all four faces of the other tet.
    // Stop whenever there are no tets left.

    for (label f = 0; f < 4; f++)
    if (T.size() > 0)
    {
        vector n = t.faceAreaNormal(f);
        n = n/Foam::mag(n);

        const scalar C = -(n & t.v()[f]);

        T = T.truncate(n,C);
    }

    return T;
}

template<class Type>
tessellation tetra::intersect(const Type& t) const
{
    return intersect(tessellation(t));
}

tessellation tetra::intersect(const tessellation& T) const
{
    tessellation R;
    forAll(T, i)
        R.append(intersect(T[i]));

    return R;
}

scalar tetra::truncationVolume(const vector& n, const scalar C) const
{
    return truncate(n,C).volume();
}

template<class Type>
scalar tetra::intersectionVolume(const Type& t) const
{
    return intersect(t).volume();
}

// Instantiate

template tessellation tetra::intersect<tetra>(const tetra&) const;
template tessellation tetra::intersect<hexa>(const hexa&) const;
template tessellation tetra::intersect<piped>(const piped&) const;
template tessellation tetra::intersect<tessellation>(const tessellation&) const;

template scalar tetra::intersectionVolume<tetra>(const tetra&) const;
template scalar tetra::intersectionVolume<hexa>(const hexa&) const;
template scalar tetra::intersectionVolume<piped>(const piped&) const;
template scalar tetra::intersectionVolume<tessellation>(const tessellation&)
    const;

}

}
