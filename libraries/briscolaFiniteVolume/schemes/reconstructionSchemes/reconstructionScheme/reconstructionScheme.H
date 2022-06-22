#ifndef reconstructionScheme_H
#define reconstructionScheme_H

#include "scheme.H"
#include "meshFields.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

template<class Type>
class reconstructionScheme
:
    public scheme
{
public:

    // Type of the reconstruction
    typedef typename outerProduct<vector,Type>::type ReconType;

    // Runtime type information
    TypeName("reconstructionScheme");

    // Declare run-time constructor selection tables
    declareRunTimeSelectionTable
    (
        autoPtr,
        reconstructionScheme,
        dictionary,
        (
            const dictionary& dict,
            const fvMesh& fvMsh
        ),
        (dict, fvMsh)
    );

    // Construct from dictionary and mesh
    reconstructionScheme(const dictionary&, const fvMesh&);

    // Copy constructor
    reconstructionScheme(const reconstructionScheme<Type>&);

    // Selector
    static autoPtr<reconstructionScheme<Type>> New(const word, const fvMesh&);

    // Destructor
    virtual ~reconstructionScheme();

    // Reconstruction of a staggered quantity into a colocated one
    virtual tmp<meshField<ReconType,colocated>> reconstruct
    (
        const meshField<Type,staggered>&
    ) = 0;
};

}

}

}

#ifdef NoRepository
    #include "reconstructionScheme.C"
#endif

#endif
