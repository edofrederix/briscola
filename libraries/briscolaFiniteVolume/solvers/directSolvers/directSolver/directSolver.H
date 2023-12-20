#ifndef directSolver_H
#define directSolver_H

#include "linearSystem.H"
#include "colocatedFields.H"
#include "staggeredFields.H"
#include "dictionary.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class SType, class Type, class MeshType>
class directSolver
{
protected:

    // Copy of the scheme dictionary
    dictionary dict_;

    // Reference to the mesh
    const fvMesh& fvMsh_;

    // Level number
    const label l_;

public:

    // Runtime type information
    TypeName("directSolver");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        directSolver,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh,
            const label l
        ),
        (dict, fvMsh, l)
    );

    // Construct from dictionary and mesh
    directSolver(const dictionary&, const fvMesh&, const label l);

    // Selector from dictionary
    static autoPtr<directSolver<SType,Type,MeshType>> New
    (
        const dictionary&,
        const fvMesh&,
        const label = 0
    );

    // Destructor
    virtual ~directSolver();

    // Access to the dictionary
    inline const dictionary& dict() const { return dict_; }

    // Access to the mesh
    inline const fvMesh& fvMsh() const { return fvMsh_; }

    // Return the mesh level
    inline label levelNum() const { return l_; }

    // Solve the given system, indicating of singular system augmentation should
    // be applied
    virtual void solve
    (
        linearSystem<SType,Type,MeshType>&,
        const List<bool>&
    ) = 0;
};

}

}

}

#ifdef NoRepository
    #include "directSolver.C"
#endif

#endif
