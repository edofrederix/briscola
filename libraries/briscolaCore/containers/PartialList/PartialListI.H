namespace Foam
{

template<class T>
inline void PartialList<T>::checkIndex(const label i) const
{
    const label j = i + startIndex_;

    if (!dataSize_)
    {
        FatalErrorInFunction
            << " attempt to access element from zero sized list"
            << abort(FatalError);
    }
    else if (j<0 || j>=dataSize_)
    {
        FatalErrorInFunction
            << "index " << j << " out of range 0 ... " << dataSize_-1
            << abort(FatalError);
    }
}

template<class T>
inline PartialList<T>::PartialList
(
    const UList<T>& list,
    const label subSize
)
:
    SubList<T>(list, subSize),
    dataSize_(list.size()),
    startIndex_(0)
{}

template<class T>
inline PartialList<T>::PartialList
(
    const UList<T>& list,
    const label subSize,
    const label startIndex
)
:
    SubList<T>(list, subSize, startIndex),
    dataSize_(list.size()),
    startIndex_(startIndex)
{}

template<class T>
inline const PartialList<T>& PartialList<T>::null()
{
    return NullObjectRef<PartialList<T>>();
}

template<class T>
inline autoPtr<PartialList<T>> PartialList<T>::clone() const
{
    return autoPtr<PartialList<T>>(new PartialList<T>(*this));
}

template<class T>
inline T& PartialList<T>::operator[](const label i)
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif

    return UList<T>::data()[i];
}

template<class T>
inline const T& PartialList<T>::operator[](const label i) const
{
    #ifdef FULLDEBUG
    checkIndex(i);
    #endif

    return UList<T>::cdata()[i];
}

}
