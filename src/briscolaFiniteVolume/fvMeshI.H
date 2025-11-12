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

            if (fvMesh::debug)
                Pout<< "Reusing " << name << endl;

            objPtr = &const_cast<Type&>(db.lookupObjectRef<Type>(name));

            // If the object is referenced by another tmp we're creating two
            // cached objects with the same name

            if (objPtr->refCount::count() > 0)
                FatalErrorInFunction
                    << "Object " << name << " is in cache but is already in use"
                    << abort(FatalError);
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
                    if (fvMesh::debug)
                        Pout<< "Reusing " << iter()->name()
                            << " for " << name << endl;

                    objPtr = const_cast<Type*>(iter());
                    objPtr->rename(name);
                }
            }

            // No available object found in cache. Allocate a new one, for which
            // we need a constructor from io and mesh.

            if (!objPtr)
            {
                if (fvMesh::debug)
                    Pout<< "Allocating " << name << endl;

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
        // Return a new unregistered object

        return
            tmp<Type>
            (
                new Type
                (
                    IOobject
                    (
                        name,
                        time().path()/"0",
                        db(),
                        IOobject::NO_READ,
                        IOobject::NO_WRITE,
                        false
                    ),
                    *this
                )
            );
    }
}

}

}

}
