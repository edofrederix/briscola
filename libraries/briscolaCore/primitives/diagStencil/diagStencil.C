#include "diagStencil.H"
#include "IOstreams.H"

const char* const Foam::diagStencil::typeName = "diagStencil";

const char* const Foam::diagStencil::componentNames[] = {"center"};

const Foam::diagStencil Foam::diagStencil::zero(diagStencil::uniform(0.0));
const Foam::diagStencil Foam::diagStencil::one(diagStencil::uniform(1.0));
const Foam::diagStencil Foam::diagStencil::max(diagStencil::uniform(vGreat));
const Foam::diagStencil Foam::diagStencil::min(diagStencil::uniform(-vGreat));
const Foam::diagStencil Foam::diagStencil::rootMax(diagStencil::uniform(rootVGreat));
const Foam::diagStencil Foam::diagStencil::rootMin(diagStencil::uniform(-rootVGreat));

Foam::word Foam::name
(
    const diagStencil& ds
)
{
    std::ostringstream buf;

    buf << '(' << name(ds.v_) << ')';

    return buf.str();
}

void Foam::writeEntry(Ostream& os, const diagStencil& value)
{
    os << value;
}

Foam::Istream& Foam::operator>>
(
    Istream& is,
    diagStencil& ds
)
{
    is.readBegin("symmStencil");

    is >> ds.v_;

    is.readEnd("symmStencil");
    is.check("operator>>(Istream&, diagStencil&)");

    return is;
}

Foam::Ostream& Foam::operator<<
(
    Ostream& os,
    const diagStencil& ds
)
{
    os << token::BEGIN_LIST << ds.v_ << token::END_LIST;

    os.check("operator<<(Ostream&, const diagStencil&)");

    return os;
}
