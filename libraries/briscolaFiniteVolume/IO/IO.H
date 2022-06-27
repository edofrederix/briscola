#ifndef IO_H
#define IO_H

#include "fvMesh.H"

#include "VectorSpace.H"
#include "CellSpace.H"

class vtkStructuredGrid;

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

    // List of stored times data (time name and number of processors)
    List<Tuple2<scalar, label>> timeData_;

    // Data file name prefix
    static word dataFileNamePrefix;

    // Name of the stored times data file
    static word timeDataFileName;

    // Default extension of the structured grid XML file
    static word XMLStructuredGridExtension;

    // Write data to given time directory for given level
    template<class MeshType>
    void writeData(const word, const label);

    // Read data from given time directory for given level
    template<class MeshType>
    void readData(const word, const label);

    // Add fields to grid for given mesh level and direction
    template<class Type, class MeshType>
    void writeFields(vtkStructuredGrid*, const label, const label) const;

    // Read fields from grid for given mesh level and direction
    template<class Type, class MeshType>
    void readFields(vtkStructuredGrid*, const label, const label) const;

    // Write

    // Write types
    template<class Type, class MeshType>
    void writeField
    (
        vtkStructuredGrid*,
        const meshDirection<Type, MeshType>&
    ) const;

    // Write scalar-like types
    template<class Type, class MeshType>
    void writeScalarField
    (
        vtkStructuredGrid*,
        const meshDirection<Type, MeshType>&
    ) const;

    // Write array-like types
    template<class Type, class MeshType>
    void writeArrayField
    (
        vtkStructuredGrid*,
        const meshDirection<Type, MeshType>&
    ) const;

    // Write two-dimensional array-like types
    template<class Type, class MeshType>
    void writeArrayArrayField
    (
        vtkStructuredGrid*,
        const meshDirection<Type, MeshType>&
    ) const;

    // Read

    // Read types
    template<class Type, class MeshType>
    void readField
    (
        vtkStructuredGrid*,
        meshDirection<Type, MeshType>&
    ) const;

    // Read scalar-like type
    template<class Type, class MeshType>
    void readScalarField
    (
        vtkStructuredGrid*,
        meshDirection<Type, MeshType>&
    ) const;

    // Write array-like types
    template<class Type, class MeshType>
    void readArrayField
    (
        vtkStructuredGrid*,
        meshDirection<Type, MeshType>&
    ) const;

    // Write two-dimensional array-like types
    template<class Type, class MeshType>
    void readArrayArrayField
    (
        vtkStructuredGrid*,
        meshDirection<Type, MeshType>&
    ) const;

    // Write VTK data to given time directory, grid and mesh level and direction
    template<class MeshType>
    void writeVTKData
    (
        const word timeName,
        vtkStructuredGrid*,
        const label,
        const label
    ) const;

    // Read VTK data for given time directory, grid and mesh level and direction
    template<class MeshType>
    void readVTKData
    (
        const word timeName,
        vtkStructuredGrid*,
        const label,
        const label
    ) const;

    // Read the stored times data for given mesh level
    void readTimeData();

    // Add given time name to time data. Replace if this time name is already
    // present
    void addToTimeData(const word);

    // Sort the time data based on time name
    void sortTimeData();

    // Write the time data to PVD files for the given level
    template<class MeshType>
    void writeTimeData(const label);

    // Return the VTK file name for given mesh level, direction and proc number
    template<class MeshType>
    word vtkFileName(const label, const label, const label) const;

    // Return the PVD file name for given mesh level and direction
    template<class MeshType>
    word pvdFileName(const label, const label) const;

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

    // Write function: checks if the current time is a write time
    template<class MeshType>
    void write(const label = 0);

    // Write now
    template<class MeshType>
    void writeNow(const label = 0);

    // Read data for current time
    template<class MeshType>
    void read(const label = 0);
};

}

}

}

#endif
