#include "arguments.H"

#include "briscolaUtils.H"

using namespace Foam;
using namespace briscola;

int main(int argc, char *argv[])
{
    arguments args(argc, argv);

    for (int i = 0; i < 6; i++)
    {
        if (i != faceNumber(faceOffsets[i]))
            FatalErrorInFunction
                << "Error in faceNumber function" << abort(FatalError);

        if (edgeNumber(faceOffsets[i]) != -1)
            FatalErrorInFunction
                << "Error in faceNumber function" << abort(FatalError);

        if (vertexNumber(faceOffsets[i]) != -1)
            FatalErrorInFunction
                << "Error in faceNumber function" << abort(FatalError);
    }

    for (int i = 0; i < 12; i++)
    {
        if (i != edgeNumber(edgeOffsets[i]))
            FatalErrorInFunction
                << "Error in edgeNumber function" << abort(FatalError);

        if (faceNumber(edgeOffsets[i]) != -1)
            FatalErrorInFunction
                << "Error in edgeNumber function" << abort(FatalError);

        if (vertexNumber(edgeOffsets[i]) != -1)
            FatalErrorInFunction
                << "Error in edgeNumber function" << abort(FatalError);
    }

    for (int i = 0; i < 8; i++)
    {
        if (i != vertexNumber(vertexOffsets[i]))
            FatalErrorInFunction
                << "Error in vertexNumber function" << abort(FatalError);

        if (faceNumber(vertexOffsets[i]) != -1)
            FatalErrorInFunction
                << "Error in vertexNumber function" << abort(FatalError);

        if (edgeNumber(vertexOffsets[i]) != -1)
            FatalErrorInFunction
                << "Error in vertexNumber function" << abort(FatalError);
    }
}
