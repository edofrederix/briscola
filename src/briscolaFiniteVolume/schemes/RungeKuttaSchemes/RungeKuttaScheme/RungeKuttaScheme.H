#ifndef RungeKuttaScheme_H
#define RungeKuttaScheme_H

#include "fvMesh.H"
#include "meshField.H"
#include "NamedEnum.H"

// Based on the work of Edo Frederix and Jannes Hopman as implemented in
// RKSymFoam, see https://github.com/janneshopman/RKSymFoam

namespace Foam
{

namespace briscola
{

namespace fv
{

class RungeKuttaScheme
:
    public regIOobject
{
public:

    // Pressure solution modes
    enum pModeType
    {
        ALLSTAGES,  // Solve for pressure every stage
        LASTSTAGE   // Solve for pressure only in the last stage
    };

    // Pressure mode names
    static const NamedEnum<pModeType,2> pModeNames;

protected:

    // Protected data

        // Reference to the mesh database
        const fvMesh& fvMsh_;

        // Copy of the Runge-Kutta dictionary
        dictionary dict_;

        // Runge-Kutta stage
        label stage_;

        // Pressure solution mode, defaulting to ALLSTAGES
        pModeType pMode_;

    // Private Member Functions

        // Disallow copy construct
        RungeKuttaScheme(const RungeKuttaScheme&);

        // Disallow default bitwise assignment
        void operator=(const RungeKuttaScheme&);

public:

    // Runtime type information
    TypeName("RungeKuttaScheme");

    // Declare run-time constructor selection table
    declareRunTimeSelectionTable
    (
        autoPtr,
        RungeKuttaScheme,
        dictionary,
        (const fvMesh& mesh),
        (mesh)
    );

    // Selector
    static autoPtr<RungeKuttaScheme> New(const fvMesh& mesh);

    // Construct from mesh
    RungeKuttaScheme(const fvMesh& fvMsh);

    // Destructor
    virtual ~RungeKuttaScheme();

    // Member Functions

        // Generate and return a list of fields to store the intermediate stage
        // sources or solutions, with given field name
        template<class Type, class MeshType>
        PtrList<meshField<Type,MeshType>> stageList(const word name) const;

        // Helper function for field types
        template<class FieldType>
        inline PtrList<FieldType> stageList(const word name) const
        {
            return
                stageList
                <
                    typename FieldType::dataType,
                    typename FieldType::meshType
                >(name);
        }

        // Access

            // Return the mesh
            inline const fvMesh& fvMsh() const;

            // Return the Runge-Kutta dictionary
            inline const dictionary& dict() const;

            // Return the pressure mode
            inline pModeType pMode() const;

            // Current stage level
            inline label stage() const;

            // Return the number of stages
            inline label numberOfStages() const;

            // Current diagonal A coefficient
            inline scalar A() const;

            // Current diagonal B coefficient
            inline scalar B() const;

            // Current C coefficient
            inline scalar C() const;

            // Return the Butcher tableau for source A
            virtual const List<scalarList>& a() const = 0;

            // Return the Butcher tableau for source B
            virtual const List<scalarList>& b() const = 0;

        // Switches

            // Explicit current stage A
            inline Switch exStageA() const;

            // Explicit current stage B
            inline Switch exStageB() const;

            // Explicit current stage
            inline Switch exStage() const;

            // Explicit scheme A
            inline Switch exA() const;

            // Explicit scheme B
            inline Switch exB() const;

            // Explicit scheme
            inline Switch ex() const;

            // Implicit current stage A
            inline Switch imStageA() const;

            // Implicit current stage B
            inline Switch imStageB() const;

            // Implicit current stage
            inline Switch imStage() const;

            // Implicit scheme A
            inline Switch imA() const;

            // Implicit scheme B
            inline Switch imB() const;

            // Explicit scheme
            inline Switch im() const;

            // First stage
            inline Switch firstStage() const;

            // Last stage
            inline Switch lastStage() const;

            // Check if stage A must be stored for future use
            inline Switch storeStageA() const;

            // Check if stage B must be stored for future use
            inline Switch storeStageB() const;

            // Solve this stage
            inline Switch solve() const;

            // Solve for pressure
            inline Switch solvePressure() const;

        // Solution control

            // Loop over stages
            inline Switch loop();

            // Sum the given list of fields for the current stage, weighted by
            // the A Butcher tableau
            template<class Type, class MeshType>
            tmp<meshField<Type,MeshType>> stageSumA
            (
                const PtrList<meshField<Type,MeshType>>& list
            ) const;

            // Sum the given list of fields for the current stage, weighted by
            // the B Butcher tableau
            template<class Type, class MeshType>
            tmp<meshField<Type,MeshType>> stageSumB
            (
                const PtrList<meshField<Type,MeshType>>& list
            ) const;

            // Sum the two given lists of fields for the current stage, weighted
            // by the A and B Butcher tableaus, respectively
            template<class Type, class MeshType>
            tmp<meshField<Type,MeshType>> stageSum
            (
                const PtrList<meshField<Type,MeshType>>& listA,
                const PtrList<meshField<Type,MeshType>>& listB
            ) const;

        // IO

            // Write data (not supported)
            bool writeData(Ostream& os) const { return os; }
};


}

}

}

#include "RungeKuttaSchemeI.H"

#endif
