#include "hexScalar.H"
#include "Istream.H"

const char* const Foam::pTraits<Foam::hexScalar>::typeName = "hexScalar";
const Foam::hexScalar Foam::pTraits<Foam::hexScalar>::zero = Foam::hexScalar::uniform(0.0);
const Foam::hexScalar Foam::pTraits<Foam::hexScalar>::one = Foam::hexScalar::uniform(1.0);
const Foam::hexScalar Foam::pTraits<Foam::hexScalar>::min = Foam::hexScalar::uniform(-vGreat);
const Foam::hexScalar Foam::pTraits<Foam::hexScalar>::max = Foam::hexScalar::uniform(vGreat);
const Foam::hexScalar Foam::pTraits<Foam::hexScalar>::rootMin = Foam::hexScalar::uniform(-rootVGreat);;
const Foam::hexScalar Foam::pTraits<Foam::hexScalar>::rootMax = Foam::hexScalar::uniform(rootVGreat);;

const char* const Foam::pTraits<Foam::hexScalar>::componentNames[] =
    { "left", "right", "bottom", "top", "aft", "fore" };

Foam::pTraits<Foam::hexScalar>::pTraits(const Foam::hexScalar& p)
:
    p_(p)
{}

Foam::pTraits<Foam::hexScalar>::pTraits(Foam::Istream& is)
{
    is >> p_;
}

template<>
const Foam::hexScalar
Foam::hexScalar::zero(hexScalar::uniform(0.0));

template<>
const Foam::hexScalar
Foam::hexScalar::one(hexScalar::uniform(1.0));

template<>
const Foam::hexScalar
Foam::hexScalar::max(hexScalar::uniform(vGreat));

template<>
const Foam::hexScalar
Foam::hexScalar::min(hexScalar::uniform(-vGreat));

template<>
const Foam::hexScalar
Foam::hexScalar::rootMax(hexScalar::uniform(rootVGreat));

template<>
const Foam::hexScalar
Foam::hexScalar::rootMin(hexScalar::uniform(-rootVGreat));
