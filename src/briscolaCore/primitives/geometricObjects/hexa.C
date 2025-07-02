#include "hexa.H"

#include "tetra.H"
#include "piped.H"
#include "tessellation.H"

namespace Foam
{

namespace briscola
{

hexa::hexa(const piped& p)
:
    v_(p.v_)
{}

scalar hexa::volume() const
{
    return tessellation(*this).volume();
}

tessellation hexa::truncate(const vector& n, const scalar C) const
{
    return tessellation(*this).truncate(n,C);
}

template<class Type>
tessellation hexa::intersect(const Type& t) const
{
    return tessellation(*this).intersect(t);
}

scalar hexa::truncationVolume(const vector& n, const scalar C) const
{
    return truncate(n,C).volume();
}

template<class Type>
scalar hexa::intersectionVolume(const Type& t) const
{
    return intersect(t).volume();
}

// Instantiate

template tessellation hexa::intersect<tetra>(const tetra&) const;
template tessellation hexa::intersect<hexa>(const hexa&) const;
template tessellation hexa::intersect<piped>(const piped&) const;
template tessellation hexa::intersect<tessellation>(const tessellation&) const;

template scalar hexa::intersectionVolume<tetra>(const tetra&) const;
template scalar hexa::intersectionVolume<hexa>(const hexa&) const;
template scalar hexa::intersectionVolume<piped>(const piped&) const;
template scalar hexa::intersectionVolume<tessellation>(const tessellation&)
    const;

}

}
