#include "stencil.H"
#include "IOstreams.H"

const char* const Foam::stencil::typeName = "stencil";

const char* const Foam::stencil::componentNames[] =
    {"center", "left", "right", "bottom", "top", "aft", "fore"};

const Foam::stencil Foam::stencil::zero(stencil::uniform(0.0));
const Foam::stencil Foam::stencil::one(stencil::uniform(1.0));
const Foam::stencil Foam::stencil::max(stencil::uniform(vGreat));
const Foam::stencil Foam::stencil::min(stencil::uniform(-vGreat));
const Foam::stencil Foam::stencil::rootMax(stencil::uniform(rootVGreat));
const Foam::stencil Foam::stencil::rootMin(stencil::uniform(-rootVGreat));

Foam::word Foam::name
(
    const stencil& st
)
{
    std::ostringstream buf;

    buf << '(' << name(st.v_[0]);

    for (direction i=1; i<7; i++)
    {
        buf << ',' << name(st.v_[i]);
    }

    buf << ')';

    return buf.str();
}

void Foam::writeEntry(Ostream& os, const stencil& value)
{
    os << value;
}

Foam::Istream& Foam::operator>>
(
    Istream& is,
    stencil& st
)
{
    is.readBegin("stencil");

    for (direction i=0; i<7; i++)
    {
        is >> st.v_[i];
    }

    is.readEnd("stencil");
    is.check("operator>>(Istream&, stencil&)");

    return is;
}

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const stencil& st
)
{
    os << token::BEGIN_LIST << st.v_[0];

    for (direction i=1; i<7; i++)
    {
        os << token::SPACE << st.v_[i];
    }

    os << token::END_LIST;
    os.check("operator<<(Ostream&, const stencil&)");

    return os;
}
