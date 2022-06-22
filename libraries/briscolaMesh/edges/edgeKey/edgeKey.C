#include "edgeKey.H"

namespace Foam
{

namespace briscola
{

edgeKey::hash::hash()
{}

edgeKey::edgeKey()
{}

edgeKey::edgeKey
(
    const label v1,
    const label v2
)
:
    Pair<label>(v1, v2)
{}

edgeKey::~edgeKey()
{}

label edgeKey::hash::operator()
(
    const edgeKey& key
) const
{
    return
        word::hash()(word(key.first()))
      + word::hash()(word(key.second()));
}

bool operator==
(
    const edgeKey& a,
    const edgeKey& b
)
{
    return mag(Pair<label>::compare(a,b));
}

bool operator!=
(
    const edgeKey& a,
    const edgeKey& b
)
{
    return !(a == b);
}

Istream& operator>>(Istream& is, edgeKey& key)
{
    const FixedList<label, 2> temp(is);

    key.first() = temp[0];
    key.second() = temp[1];

    return is;
}

Ostream& operator<<(Ostream& os, const edgeKey& key)
{
    os  << token::BEGIN_LIST
        << key.first()
        << token::SPACE
        << key.second()
        << token::END_LIST;

    return os;
}

}

}
