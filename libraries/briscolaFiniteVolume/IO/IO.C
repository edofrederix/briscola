#include "IO.H"
#include "OSspecific.H"
#include "Time.H"
#include "SortableList.H"
#include "colocatedFields.H"
#include "staggeredFields.H"
#include "floatVector.H"

namespace Foam
{

namespace briscola
{

namespace fv
{

defineTypeNameAndDebug(IO, 0);

word IO::dataFileNamePrefix("briscola");
word IO::vtkHeader("# vtk DataFile Version 3.0");

template<class MeshType>
void IO::writeData(const word timeName, const label l)
{
    // Prepare points

    const meshLevel<vector,MeshType>& ccl =
        fvMsh_.template metrics<MeshType>().cellCenters()[l];

    const partLevelPoints& points = fvMsh_[l].points();

    const bool ascii = (fvMsh_.time().writeFormat() == IOstream::ASCII);

    forAll(ccl, d)
    {
        const word vtkFile = vtkFileName<MeshType>(l,d);

        const fileName vtkFilePath
        (
            fvMsh_.time().path()/timeName/vtkFile
        );

        // Create file stream on master

        autoPtr<std::ofstream> filePtr;

        if (Pstream::master())
            filePtr.reset(new std::ofstream(vtkFilePath));

        // Write header

        if (Pstream::master())
        {
            std::ofstream& file = filePtr();

            file<< vtkHeader << std::endl
                << fvMsh_.time().caseName() << " "
                << Pstream::nProcs() << std::endl;

            if (ascii)
            {
                file<< "ASCII" << std::endl;
            }
            else
            {
                file<< "BINARY" << std::endl;
            }

            file<< "DATASET UNSTRUCTURED_GRID" << std::endl;
        }

        // Create points data

        const labelVector S(ccl[d].S());
        const labelVector E(ccl[d].E());

        List<floatScalar> myPoints(cmptProduct(E-S+unitXYZ)*3);

        const vector shift(MeshType::shift[d]);

        label c = 0;
        for (int i = S.x(); i < E.x()+1; i++)
        for (int j = S.y(); j < E.y()+1; j++)
        for (int k = S.z(); k < E.z()+1; k++)
        {
            const floatVector p =
                floatVector(points.interp(vector(i,j,k)+shift));

            myPoints[c++] = p.x();
            myPoints[c++] = p.y();
            myPoints[c++] = p.z();
        }

        // Write points data

        List<labelVector> shapes(Pstream::nProcs());
        shapes[Pstream::myProcNo()] = E-S;
        Pstream::gatherList(shapes);

        if (Pstream::master())
        {
            label nPoints = 0;

            forAll(shapes, proc)
                nPoints += cmptProduct(shapes[proc]+unitXYZ);

            std::ofstream& file = filePtr();

            file<< "POINTS " << nPoints << " float" << std::endl;
        }

        IO::writeList
        (
            filePtr,
            myPoints,
            word("points"),
            1,
            fvMsh_.time().writeFormat() == IOstream::ASCII,
            l*MeshType::numberOfDirections + d,
            false
        );

        // Write cell data

        if (Pstream::master())
        {
            label nCells = 0;

            forAll(shapes, proc)
                nCells += cmptProduct(shapes[proc]);

            std::ofstream& file = filePtr();

            file<< "CELLS " << nCells << " " << 9*nCells << std::endl;

            label cursor = 0;

            forAll(shapes, proc)
            {
                const label l = shapes[proc].x();
                const label m = shapes[proc].y();
                const label n = shapes[proc].z();

                const label L = l+1;
                const label M = m+1;
                const label N = n+1;

                labelList buffer(l*m*n*9);

                label c = 0;

                for(int i = 0; i < l; i++)
                for(int j = 0; j < m; j++)
                for(int k = 0; k < n; k++)
                {
                    buffer[c++] = 8;
                    buffer[c++] = (cursor + (i  )*M*N + (j  )*N+k);
                    buffer[c++] = (cursor + (i+1)*M*N + (j  )*N+k);
                    buffer[c++] = (cursor + (i+1)*M*N + (j+1)*N+k);
                    buffer[c++] = (cursor + (i  )*M*N + (j+1)*N+k);
                    buffer[c++] = (cursor + (i  )*M*N + (j  )*N+k+1);
                    buffer[c++] = (cursor + (i+1)*M*N + (j  )*N+k+1);
                    buffer[c++] = (cursor + (i+1)*M*N + (j+1)*N+k+1);
                    buffer[c++] = (cursor + (i  )*M*N + (j+1)*N+k+1);
                }

                if (ascii)
                {
                    forAll(buffer, i)
                    {
                        file<< buffer[i] << " ";
                    }
                }
                else
                {
                    #ifdef LITTLEENDIAN
                    swapWords
                    (
                        buffer.size(),
                        reinterpret_cast<label*>(buffer.begin())
                    );
                    #endif

                    file.write
                    (
                        reinterpret_cast<char*>(buffer.begin()),
                        buffer.byteSize()
                    );
                }

                cursor += L*M*N;
            }

            file<< std::endl;

            file<< "CELL_TYPES " << nCells << std::endl;

            if (ascii)
            {
                for (int i = 0; i < nCells; i++)
                {
                    file<< "12 ";
                }
            }
            else
            {
                labelList buffer(nCells, 12);

                #ifdef LITTLEENDIAN
                swapWords
                (
                    buffer.size(),
                    reinterpret_cast<label*>(buffer.begin())
                );
                #endif

                file.write
                (
                    reinterpret_cast<char*>(buffer.begin()),
                    buffer.byteSize()
                );
            }

            file<< std::endl;
            file<< "CELL_DATA " << nCells << std::endl;
        }

        // Write fields as cell data. Order is important.

        writeFields<label,MeshType>(filePtr, l, d);
        writeFields<scalar,MeshType>(filePtr, l, d);

        writeFields<vector,MeshType>(filePtr, l, d);
        writeFields<tensor,MeshType>(filePtr, l, d);
        writeFields<diagTensor,MeshType>(filePtr, l, d);
        writeFields<sphericalTensor,MeshType>(filePtr, l, d);
        writeFields<symmTensor,MeshType>(filePtr, l, d);

        writeFields<faceScalar,MeshType>(filePtr, l, d);
        writeFields<faceVector,MeshType>(filePtr, l, d);
    }
}

template<class MeshType>
void IO::readData(const word timeName, const label l)
{
    const meshLevel<vector,MeshType>& ccl =
        fvMsh_.template metrics<MeshType>().cellCenters()[l];

    for (label d = 0; d < MeshType::numberOfDirections; d++)
    {
        const word vtkFile = vtkFileName<MeshType>(l,d);

        const fileName vtkFilePath
        (
            fvMsh_.time().path()/timeName/vtkFile
        );

        // Create file stream

        autoPtr<std::ifstream> filePtr;

        filePtr.reset(new std::ifstream(vtkFilePath));

        std::ifstream& file = filePtr();

        if (file.fail())
        {
            return;
        }

        // Read header

        char buffer[256];

        file.getline(buffer, sizeof(buffer));

        word header(buffer);

        if (header != vtkHeader && Pstream::master())
        {
            FatalErrorInFunction
                << "Invalid header in file " << vtkFilePath
                << endl << abort(FatalError);
        }

        // Case name and number of procs

        word dummy;
        label nProcs;

        file>> dummy >> nProcs;
        nextLine(file);

        if (nProcs != Pstream::nProcs() && Pstream::master())
            FatalErrorInFunction
                << "The vtk file " << vtkFilePath << nl
                << "was generated by " << nProcs << " processors "
                << "which is not compatible with the current number "
                << "of processors (" << Pstream::nProcs() << ")"
                << endl << abort(FatalError);

        // Ascii or binary

        word dataType;

        file>> dataType;
        nextLine(file);

        bool ascii = false;

        if (dataType == "ASCII")
        {
            ascii = true;
        }
        else if (dataType == "BINARY")
        {
            ascii = false;
        }
        else if (Pstream::master())
        {
            FatalErrorInFunction
                << "Could not determine if file " << vtkFilePath << nl
                << "is binary or ascii"
                << endl << abort(FatalError);
        }

        // Skip dataset

        nextLine(file);

        // Skip points

        if (ascii)
        {
            nextLine(file);
            nextLine(file);
        }
        else
        {
            word dummy;
            label nPoints;

            file>> dummy >> nPoints;

            nextLine(file);
            file.ignore(nPoints*3*sizeof(floatScalar));
            nextLine(file);
        }

        // Skip cells

        if (ascii)
        {
            nextLine(file);
            nextLine(file);
        }
        else
        {
            word dummy;
            label n;

            file>> dummy >> dummy >> n;

            nextLine(file);
            file.ignore(n*sizeof(label));
            nextLine(file);
        }

        // Skip cell types

        if (ascii)
        {
            nextLine(file);
            nextLine(file);
        }
        else
        {
            word dummy;
            label nCells;

            file>> dummy >> nCells;

            nextLine(file);
            file.ignore(nCells*sizeof(label));
            nextLine(file);
        }

        // Read number of cells

        label nCells;

        file>> dummy >> nCells;

        labelList sizes(Pstream::nProcs());
        sizes[Pstream::myProcNo()] = ccl[d].size();
        Pstream::gatherList(sizes);
        Pstream::scatterList(sizes);

        if (sum(sizes) != nCells && Pstream::master())
            FatalErrorInFunction
                << "Trying to read " << sum(sizes) << " cells "
                << "but found data for " << nCells << " cells in " << nl
                << vtkFilePath << endl << abort(FatalError);

        // Read fields from cell data. Order is important.

        readFields<label,MeshType>(filePtr, ascii, l, d);
        readFields<scalar,MeshType>(filePtr, ascii, l, d);

        readFields<vector,MeshType>(filePtr, ascii, l, d);
        readFields<tensor,MeshType>(filePtr, ascii, l, d);
        readFields<diagTensor,MeshType>(filePtr, ascii, l, d);
        readFields<sphericalTensor,MeshType>(filePtr, ascii, l, d);
        readFields<symmTensor,MeshType>(filePtr, ascii, l, d);

        readFields<faceScalar,MeshType>(filePtr, ascii, l, d);
        readFields<faceVector,MeshType>(filePtr, ascii, l, d);
    }
}

template<class Type, class MeshType>
void IO::writeFields
(
    autoPtr<std::ofstream>& filePtr,
    const label l,
    const label d
) const
{
    typedef meshField<Type,MeshType> FieldType;

    const wordList names(fvMsh_.time().names<FieldType>());

    label count = 0;

    forAll(names, i)
    {
        const FieldType& field =
            fvMsh_.time().lookupObject<FieldType>(names[i]);

        if (field.writeOpt() == IOobject::AUTO_WRITE)
        {
            count++;
        }
    }

    if (Pstream::master())
    {
        std::ofstream& file = filePtr();

        file<< "FIELD " << FieldType::typeName << "s "
            << count << std::endl;
    }

    if (count > 0)
    {
        forAll(names, i)
        {
            const FieldType& field =
                fvMsh_.time().lookupObject<FieldType>(names[i]);

            if (field.writeOpt() == IOobject::AUTO_WRITE)
            {
                writeField(filePtr, field[l][d]);
            }
        }
    }
}

template<class Type, class MeshType>
void IO::readFields
(
    autoPtr<std::ifstream>& filePtr,
    const bool ascii,
    const label l,
    const label d
) const
{
    typedef meshField<Type, MeshType> FieldType;

    word dummy, type;
    label nf;

    std::ifstream& file = filePtr();

    file>> dummy >> type >> nf;

    if (word(FieldType::typeName)+"s" != type && Pstream::master())
    {
        FatalErrorInFunction
            << type << endl
            << "Could not find data for field type " << FieldType::typeName
            << endl << abort(FatalError);
    }

    if (nf > 0)
    {
        word name;

        for (int i = 0; i < nf; i++)
        {
            label nComponents;
            label nCells;

            file>> name >> nComponents >> nCells;
            nextLine(file);

            if (fvMsh_.time().foundObject<FieldType>(name))
            {
                FieldType& field =
                    fvMsh_.time().lookupObjectRef<FieldType>(name);

                if
                (
                    field.readOpt() == IOobject::MUST_READ
                 || field.readOpt() == IOobject::READ_IF_PRESENT
                )
                {
                    readField(filePtr, ascii, field[l][d]);
                }
                else
                {
                    if (ascii)
                    {
                        nextLine(file);
                    }
                    else
                    {
                        file.ignore(nCells*nComponents*sizeof(floatScalar));
                    }
                }
            }
            else
            {
                // Ignore. Alternatively we could create a field and store in
                // the registry.

                if (ascii)
                {
                    nextLine(file);
                }
                else
                {
                    file.ignore(nCells*nComponents*sizeof(floatScalar));
                }
            }
        }
    }
}

template<class MeshType>
word IO::vtkFileName
(
    const label l,
    const label d
) const
{
    return
        dataFileNamePrefix
      + "_"
      + (l > 0 ? word(Foam::name(l) + "_") : "")
      + MeshType::typeName
      + (
            MeshType::numberOfDirections > 1
          ? word("_" + word(MeshType::directionNames[d]))
          : word("")
        )
      + word(".vtk");
}

template<class MeshType>
word IO::seriesFileName(const label l, const label d) const
{
    return
        dataFileNamePrefix
      + "_"
      + (l > 0 ? word(Foam::name(l) + "_") : "")
      + MeshType::typeName
      + (
            MeshType::numberOfDirections > 1
          ? word("_" + word(MeshType::directionNames[d]))
          : word("")
        )
      + ".vtk.series";
}

template<class MeshType>
void IO::writeSeries(const label l) const
{
    instantList times(fvMsh_.time().times());

    for (label d = 0; d < MeshType::numberOfDirections; d++)
    {
        OFstream seriesFile(fvMsh_.time().path()/seriesFileName<MeshType>(l,d));

        seriesFile
            << "{" << nl
            << "    \"file-series-version\" : \"1.0\"," << nl
            << "    \"files\" :" << nl
            << "    [" << nl;

        forAll(times, i)
        if (times[i].value() != 0)
        {
            seriesFile
                << "        { \"name\" : \""
                << times[i].name() << "/" << vtkFileName<MeshType>(l,d)
                << "\", \"time\" : "
                << times[i].name() << " }," << nl;
        }

        seriesFile
            << "    ]" << nl
            << "}" << endl;
    }
}

void IO::swapWord(label& word32) const
{
    char* mem = reinterpret_cast<char*>(&word32);

    char a = mem[0];
    mem[0] = mem[3];
    mem[3] = a;

    a = mem[1];
    mem[1] = mem[2];
    mem[2] = a;
}


void IO::swapWords(const label nWords, label* words32) const
{
    for (label i = 0; i < nWords; i++)
    {
        swapWord(words32[i]);
    }
}

IO::IO(const fvMesh& fvMsh)
:
    fvMsh_(fvMsh)
{}

IO::~IO()
{}

template<class MeshType>
void IO::write(const label l)
{
    if (fvMsh_.time().writeTime())
    {
        writeNow<MeshType>(l);
    }
}

template<class MeshType>
void IO::writeNow(const label l)
{
    const Time& time = fvMsh_.time();

    const word timeName(time.timeName());

    if (Pstream::master())
    {
        if (!isDir(time.path()/timeName))
        {
            mkDir(time.path()/timeName);
        }
    }

    // Wait

    if (Pstream::parRun())
        returnReduce(1.0,sumOp<scalar>());

    // Write data

    writeData<MeshType>(timeName,l);

    // Write series

    writeSeries<MeshType>(l);
}

template<class MeshType>
void IO::read(const label l)
{
    const Time& time = fvMsh_.time();

    const word timeName(time.timeName());

    if (Pstream::master())
    {
        if (!isDir(time.path()/timeName) || timeName == "0")
        {
            FatalErrorInFunction
                << "Cannot read from time " << timeName << endl
                << exit(FatalError);
        }
    }

    // Read data

    readData<MeshType>(timeName,l);
}

// Instantiate

template void IO::read<colocated>(const label);
template void IO::read<staggered>(const label);

template void IO::write<colocated>(const label);
template void IO::write<staggered>(const label);

}

}

}
