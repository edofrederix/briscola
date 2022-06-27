#include "CellSpace.H"
#include "IOstreams.H"

#include <sstream>

namespace Foam
{

template<class Form, class Cmpt, direction Ncmpts>
CellSpace<Form, Cmpt, Ncmpts>::CellSpace
(
    Istream& is
)
{
    // Read beginning of CellSpace<Cmpt>
    is.readBegin("CellSpace<Form, Cmpt, Ncmpts>");

    for (direction i=0; i<Ncmpts; i++)
    {
        is >> v_[i];
    }

    // Read end of CellSpace<Cmpt>
    is.readEnd("CellSpace<Form, Cmpt, Ncmpts>");

    // Check state of Istream
    is.check("CellSpace<Form, Cmpt, Ncmpts>::CellSpace(Istream&)");
}

template<class Form, class Cmpt, direction Ncmpts>
word name
(
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    std::ostringstream buf;

    buf << '(' << cs.v_[0];

    for (direction i=1; i<Ncmpts; i++)
    {
        buf << ',' << cs.v_[i];
    }

    buf << ')';

    return buf.str();
}

template<class Form, class Cmpt, direction Ncmpts>
void writeEntry(Ostream& os, const CellSpace<Form, Cmpt, Ncmpts>& value)
{
    os << value;
}

template<class Form, class Cmpt, direction Ncmpts>
Istream& operator>>
(
    Istream& is,
    CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    // Read beginning of CellSpace<Form, Cmpt, Ncmpts>
    is.readBegin("CellSpace<Form, Cmpt, Ncmpts>");

    for (direction i=0; i<Ncmpts; i++)
    {
        is >> cs.v_[i];
    }

    // Read end of CellSpace<Form, Cmpt, Ncmpts>
    is.readEnd("CellSpace<Form, Cmpt, Ncmpts>");

    // Check state of Istream
    is.check("operator>>(Istream&, CellSpace<Form, Cmpt, Ncmpts>&)");

    return is;
}

template<class Form, class Cmpt, direction Ncmpts>
Ostream& operator<<
(
    Ostream& os,
    const CellSpace<Form, Cmpt, Ncmpts>& cs
)
{
    os << token::BEGIN_LIST << cs.v_[0];

    for (direction i=1; i<Ncmpts; i++)
    {
        os << token::SPACE << cs.v_[i];
    }

    os << token::END_LIST;

    // Check state of Ostream
    os.check("operator<<(Ostream&, const CellSpace<Form, Cmpt, Ncmpts>&)");

    return os;
}

}