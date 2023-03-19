#ifndef IO_H
#define IO_H

#include "fvMesh.H"

#include "IFstream.H"
#include "OFstream.H"
#include "VectorSpace.H"
#include "CellSpace.H"

#if defined(LITTLE_ENDIAN) \
 || defined(_LITTLE_ENDIAN) \
 || defined(__LITTLE_ENDIAN)
    #define LITTLEENDIAN 1
#elif defined(BIG_ENDIAN) || defined(_BIG_ENDIAN) || defined(__BIG_ENDIAN)
    #undef LITTLEENDIAN
#else
    #error "Cannot find LITTLE_ENDIAN or BIG_ENDIAN symbol defined."
    #error "Please add to compilation options"
#endif

namespace Foam
{

namespace briscola
{

namespace fv
{

class fvMesh;

class IO
{
    // Reference to mesh
    const fvMesh& fvMsh_;

    // Write processor data in partitioned mode. Can be set in controlDict with
    // the 'writePartitioned' option. Defaults to false.
    bool partitioned_;

    // Write ghosts. Can be set in controlDict with the writeGhosts option.
    // Defaults to false.
    bool ghosts_;

    // Data file name prefix
    static word dataFileNamePrefix;

    // VTK header
    static word vtkHeader;

    // Write data to given time directory for given level
    template<class MeshType>
    void writeData(const word, const label);

    // Read data from given time directory for given level
    template<class MeshType>
    void readData(const word, const label);

    // Write fields of given type, level and direction
    template<class Type, class MeshType>
    void writeFields
    (
        List<autoPtr<std::ofstream>>&,
        const label,
        const label
    ) const;

    // Read fields of given type, level and direction
    template<class Type, class MeshType>
    void readFields
    (
        autoPtr<std::ifstream>&,
        const bool,
        const label,
        const label
    ) const;

    // Write

    // Basic float data list writer for a file stream, list of floats, field
    // name, number of components, ascii bool, tag number and write header bool.
    void writeList
    (
        List<autoPtr<std::ofstream>>&,
        List<floatScalar>&,
        const word,
        const label,
        const bool,
        const label,
        const bool = true
    ) const;

    // Write types
    template<class Type, class MeshType>
    void writeField
    (
        List<autoPtr<std::ofstream>>&,
        const meshDirection<Type, MeshType>&
    ) const;

    // Write scalar-like types
    template<class Type, class MeshType>
    void writeScalarField
    (
        List<autoPtr<std::ofstream>>&,
        const meshDirection<Type, MeshType>&
    ) const;

    // Write array-like types
    template<class Type, class MeshType>
    void writeArrayField
    (
        List<autoPtr<std::ofstream>>&,
        const meshDirection<Type, MeshType>&
    ) const;

    // Write two-dimensional array-like types
    template<class Type, class MeshType>
    void writeArrayArrayField
    (
        List<autoPtr<std::ofstream>>&,
        const meshDirection<Type, MeshType>&
    ) const;

    // Read

    // Basic float data list reader for a file stream, list of floats, ascii
    // bool and tag number
    void readList
    (
        autoPtr<std::ifstream>&,
        List<floatScalar>&,
        const bool,
        const label
    ) const;

    // Read types
    template<class Type, class MeshType>
    void readField
    (
        autoPtr<std::ifstream>&,
        const bool,
        meshDirection<Type, MeshType>&
    ) const;

    // Read scalar-like type
    template<class Type, class MeshType>
    void readScalarField
    (
        autoPtr<std::ifstream>&,
        const bool,
        meshDirection<Type, MeshType>&
    ) const;

    // Write array-like types
    template<class Type, class MeshType>
    void readArrayField
    (
        autoPtr<std::ifstream>&,
        const bool,
        meshDirection<Type, MeshType>&
    ) const;

    // Write two-dimensional array-like types
    template<class Type, class MeshType>
    void readArrayArrayField
    (
        autoPtr<std::ifstream>&,
        const bool,
        meshDirection<Type, MeshType>&
    ) const;

    // Return the VTK file name for given mesh level, direction and proc
    template<class MeshType>
    word vtkFileName(const label, const label, const label) const;

    // Return the VTK series file name for given mesh level, direction and proc
    template<class MeshType>
    word seriesFileName(const label, const label, const label) const;

    // Write the series file for given mesh level
    template<class MeshType>
    void writeSeries(const label) const;

    // Swap halves of word (copied from vtkWriteOps.H)
    void swapWord(label& word32) const;

    // Swap halves of word (copied from vtkWriteOps.H)
    void swapWords(const label nWords, label* words32) const;

    // Jump to the next line in an input file stream
    inline void nextLine(std::ifstream& file) const
    {
        file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    }

public:

    // Declare name of the class and its debug switch
    ClassName("IO");

    // Construct from mesh
    IO(const fvMesh&);

    // Disallow default bitwise copy construction
    IO(const IO&);

    // Disallow default bitwise assignment
    void operator=(const IO&);

    // Destructor
    ~IO();

    // Return the mesh
    inline const fvMesh& fvMsh() const
    {
        return fvMsh_;
    }

    // Return access to the partitioned io bool
    inline bool& partitioned() { return partitioned_; }
    inline bool partitioned() const { return partitioned_; }

    // Return access to the ghosts io bool
    inline bool& ghosts() { return ghosts_; }
    inline bool ghosts() const { return ghosts_; }

    // Enable partitioned io
    inline void enablePartitioned()
    {
        partitioned_ = true;
    }

    // Disable partitioned io
    inline void disablePartitioned()
    {
        partitioned_ = false;
    }

    // Enable ghost cell io
    inline void enableGhosts()
    {
        ghosts_ = true;
    }

    // Disable ghost cell io
    inline void disableGhosts()
    {
        ghosts_ = false;
    }

    // Write function for given mesh level
    template<class MeshType>
    void write(const label = 0);

    // Write now for given mesh level
    template<class MeshType>
    void writeNow(const label = 0);

    // Read data for current time and given mesh level
    template<class MeshType>
    void read(const label = 0);
};

}

}

}

#endif
