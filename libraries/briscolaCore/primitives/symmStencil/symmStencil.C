#include "symmStencil.H"
#include "IOstreams.H"

const char* const Foam::symmStencil::typeName = "symmStencil";

const char* const Foam::symmStencil::componentNames[] = {"center", "left", "bottom", "aft"};

const Foam::symmStencil Foam::symmStencil::zero(symmStencil::uniform(0.0));
const Foam::symmStencil Foam::symmStencil::one(symmStencil::uniform(1.0));
const Foam::symmStencil Foam::symmStencil::max(symmStencil::uniform(vGreat));
const Foam::symmStencil Foam::symmStencil::min(symmStencil::uniform(-vGreat));
const Foam::symmStencil Foam::symmStencil::rootMax(symmStencil::uniform(rootVGreat));
const Foam::symmStencil Foam::symmStencil::rootMin(symmStencil::uniform(-rootVGreat));

Foam::word Foam::name
(
    const symmStencil& ss
)
{
    std::ostringstream buf;

    buf << '(' << name(ss.v_[0]);

    for (direction i=1; i<4; i++)
    {
        buf << ',' << name(ss.v_[i]);
    }

    buf << ')';

    return buf.str();
}

void Foam::writeEntry(Ostream& os, const symmStencil& value)
{
    os << value;
}

Foam::Istream& Foam::operator>>
(
    Istream& is,
    symmStencil& ss
)
{
    is.readBegin("symmStencil");

    for (direction i=0; i<4; i++)
    {
        is >> ss.v_[i];
    }

    is.readEnd("symmStencil");
    is.check("operator>>(Istream&, symmStencil&)");

    return is;
}

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const symmStencil& ss
)
{
    os << token::BEGIN_LIST << ss.v_[0];

    for (direction i=1; i<4; i++)
    {
        os << token::SPACE << ss.v_[i];
    }

    os << token::END_LIST;
    os.check("operator<<(Ostream&, const symmStencil&)");

    return os;
}
