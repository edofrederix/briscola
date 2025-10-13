namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
tmp<Type> fvMesh::NewObject(const word name, const bool cached) const
{
    if (cached)
    {
        // Cached objects are stored in the special cache registry

        const objectRegistry& db = cache();

        Type* objPtr = nullptr;

        if (db.foundObject<Type>(name))
        {
            // The object was found in the cache

            Info<< "Reusing " << name << endl;

            objPtr = &const_cast<Type&>(db.lookupObjectRef<Type>(name));

            // If the object is referenced by another tmp we're creating two
            // cached objects with the same name

            if (objPtr->refCount::count() > 0)
                FatalErrorInFunction
                    << "Object " << name << " is in cache but is already in use"
                    << abort(FatalError);

            #ifndef NO_BLOCK_ZERO_INIT
            *objPtr = Zero;
            #endif
        }
        else
        {
            // Look in the registry for another available field

            HashTable<const Type*> objects(db.lookupClass<Type>());

            forAllIter
            (
                typename HashTable<const Type*>,
                objects,
                iter
            )
            {
                if (iter()->refCount::count() < 1)
                {
                    Info<< "Reusing " << iter()->name()
                        << " for " << name << endl;

                    objPtr = const_cast<Type*>(iter());
                    objPtr->rename(name);

                    #ifndef NO_BLOCK_ZERO_INIT
                    *objPtr = Zero;
                    #endif
                }
            }

            // No available object found in cache. Allocate a new one, for which
            // we need a constructor from io and mesh.

            if (!objPtr)
            {
                Info<< "Allocating " << name << endl;

                objPtr =
                    new Type
                    (
                        IOobject
                        (
                            name,
                            time().path()/"0",
                            cache()
                        ),
                        *this
                    );

                // Store in registry
                db.store(objPtr);

                // Set object as cached
                objPtr->cached() = true;
            }
        }

        // Create tmp from pointer
        tmp<Type> tObj(objPtr);

        // Increase reference count for the tmp
        objPtr->operator++();

        return tObj;
    }
    else
    {
        return
            tmp<Type>
            (
                new Type
                (
                    IOobject
                    (
                        name,
                        time().path()/"0",
                        db()
                    ),
                    *this
                )
            );
    }
}

}

}

}
