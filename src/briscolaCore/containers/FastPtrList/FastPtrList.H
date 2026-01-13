#ifndef FastPtrList_H
#define FastPtrList_H

#include "PtrList.H"

namespace Foam
{

template<class T>
class FastPtrList
:
    public PtrList<T>
{
public:

    // Use PtrList constructors
    using PtrList<T>::PtrList;

    // Override operator[] to avoid valid pointer check

    inline T& operator[](const label i)
    {
        #ifdef FULLDEBUG

        if (i < 0 || i >= PtrList<T>::size())
            FatalErrorInFunction
                << "Index " << i << " out of range 0 ... "
                << PtrList<T>::size()-1
                << endl << abort(FatalError);

        if (!PtrList<T>::operator()(i))
        {
            FatalErrorInFunction
                << "hanging pointer at index " << i
                << " (size " << PtrList<T>::size()
                << "), cannot dereference"
                << abort(FatalError);
        }

        #endif

        return *const_cast<T*>(PtrList<T>::operator()(i));
    }

    inline const T& operator[](const label i) const
    {
        #ifdef FULLDEBUG

        if (i < 0 || i >= PtrList<T>::size())
            FatalErrorInFunction
                << "Index " << i << " out of range 0 ... "
                << PtrList<T>::size()-1
                << endl << abort(FatalError);

        if (!PtrList<T>::operator()(i))
        {
            FatalErrorInFunction
                << "hanging pointer at index " << i
                << " (size " << PtrList<T>::size()
                << "), cannot dereference"
                << abort(FatalError);
        }

        #endif

        return *PtrList<T>::operator()(i);
    }
};

}

#endif
