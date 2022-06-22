#include "hexVector.H"
#include "Istream.H"

const char* const Foam::pTraits<Foam::hexVector>::typeName = "hexVector";
const Foam::hexVector Foam::pTraits<Foam::hexVector>::zero = Foam::hexVector::uniform(vector::zero);
const Foam::hexVector Foam::pTraits<Foam::hexVector>::one = Foam::hexVector::uniform(vector::one);
const Foam::hexVector Foam::pTraits<Foam::hexVector>::min = Foam::hexVector::uniform(-vector::max);
const Foam::hexVector Foam::pTraits<Foam::hexVector>::max = Foam::hexVector::uniform(vector::min);
const Foam::hexVector Foam::pTraits<Foam::hexVector>::rootMin = Foam::hexVector::uniform(-vector::rootMax);;
const Foam::hexVector Foam::pTraits<Foam::hexVector>::rootMax = Foam::hexVector::uniform(vector::rootMin);;

const char* const Foam::pTraits<Foam::hexVector>::componentNames[] =
    { "left", "right", "bottom", "top", "aft", "fore" };

Foam::pTraits<Foam::hexVector>::pTraits(const Foam::hexVector& p)
:
    p_(p)
{}

Foam::pTraits<Foam::hexVector>::pTraits(Foam::Istream& is)
{
    is >> p_;
}

template<>
const Foam::hexVector
Foam::hexVector::zero(hexVector::uniform(vector::zero));

template<>
const Foam::hexVector
Foam::hexVector::one(hexVector::uniform(vector::one));

template<>
const Foam::hexVector
Foam::hexVector::max(hexVector::uniform(vector::max));

template<>
const Foam::hexVector
Foam::hexVector::min(hexVector::uniform(vector::min));

template<>
const Foam::hexVector
Foam::hexVector::rootMax(hexVector::uniform(vector::rootMax));

template<>
const Foam::hexVector
Foam::hexVector::rootMin(hexVector::uniform(vector::rootMin));
