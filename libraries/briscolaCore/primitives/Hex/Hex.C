#include "Hex.H"

#include "IOstreams.H"
#include <sstream>

template<class Type>
const char* const Foam::Hex<Type>::componentNames[] =
    {"left", "right", "bottom", "top", "aft", "fore"};

template<class Type>
Foam::word Foam::name
(
    const Hex<Type>& h
)
{
    std::ostringstream buf;

    buf << '(' << name(h.v_[0]);

    for (direction i=1; i<3; i++)
    {
        buf << ',' << name(h.v_[i]);
    }

    buf << ')';

    return buf.str();
}

template<class Type>
void Foam::writeEntry(Ostream& os, const Hex<Type>& h)
{
    os << h;
}

template<class Type>
Foam::Istream& Foam::operator>>
(
    Istream& is,
    Hex<Type>& h
)
{
    is.readBegin("Hex");

    for (direction i=0; i < Hex<Type>::size(); i++)
    {
        is >> h.v_[i];
    }

    is.readEnd("Hex");
    is.check("operator>>(Istream&, Hex&)");

    return is;
}

template<class Type>
Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const Hex<Type>& h
)
{
    os << token::BEGIN_LIST << h.v_[0];

    for (direction i=1; i < Hex<Type>::size(); i++)
    {
        os << token::SPACE << h.v_[i];
    }

    os << token::END_LIST;
    os.check("operator<<(Ostream&, const Hex&)");

    return os;
}
