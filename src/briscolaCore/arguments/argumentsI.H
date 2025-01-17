#include "arguments.H"

namespace Foam
{

namespace briscola
{

// Template specialization for string
template<>
inline string
arguments::argRead<string>(const label index) const
{
    return args_[index];
}

// Template specialization for word
template<>
inline word
arguments::argRead<word>(const label index) const
{
    return args_[index];
}

// Template specialization for fileName
template<>
inline fileName
arguments::argRead<fileName>(const label index) const
{
    return args_[index];
}

// Template specialization for string
template<>
inline string
arguments::optionRead<string>(const word& opt) const
{
    return options_[opt];
}

// Template specialization for word
template<>
inline word
arguments::optionRead<word>(const word& opt) const
{
    return options_[opt];
}

// Template specialization for fileName
template<>
inline fileName
arguments::optionRead<fileName>(const word& opt) const
{
    return options_[opt];
}

template<class T>
inline T arguments::argRead(const label index) const
{
    T val;

    IStringStream(args_[index])() >> val;
    return val;
}


template<class T>
inline T arguments::optionRead(const word& opt) const
{
    T val;

    optionLookup(opt)() >> val;
    return val;
}


template<class T>
inline bool arguments::optionReadIfPresent
(
    const word& opt,
    T& val
) const
{
    if (optionFound(opt))
    {
        val = optionRead<T>(opt);
        return true;
    }
    else
    {
        return false;
    }
}

template<class T>
inline bool arguments::optionReadIfPresent
(
    const word& opt,
    T& val,
    const T& deflt
) const
{
    if (optionReadIfPresent<T>(opt, val))
    {
        return true;
    }
    else
    {
        val = deflt;
        return false;
    }
}

template<class T>
inline T arguments::optionLookupOrDefault
(
    const word& opt,
    const T& deflt
) const
{
    if (optionFound(opt))
    {
        return optionRead<T>(opt);
    }
    else
    {
        return deflt;
    }
}

inline const string& arguments::operator[](const label index) const
{
    return args_[index];
}

inline const string& arguments::operator[](const word& opt) const
{
    return options_[opt];
}

}

}
