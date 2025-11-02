#ifndef RungeKuttaScheme_H
#define RungeKuttaScheme_H

#include "fvMesh.H"
#include "meshField.H"
#include "FastPtrList.H"
#include "scalarMatrices.H"

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
protected:

    // Protected data

        // Reference to the mesh database
        const fvMesh& fvMsh_;

        // Copy of the Runge-Kutta dictionary
        dictionary dict_;

        // Runge-Kutta stage
        label stage_;

        // Number of stages
        label nStages_;

        // Number of time steps
        label nSteps_;

        // Butcher tableaus
        scalarRectangularMatrix a_;
        scalarRectangularMatrix b_;

        // Stage sources
        FastPtrList<FastPtrList<colocatedScalarField>> coloScalarStageSources_;
        FastPtrList<FastPtrList<colocatedVectorField>> coloVectorStageSources_;
        FastPtrList<FastPtrList<staggeredScalarField>> stagScalarStageSources_;
        FastPtrList<FastPtrList<staggeredVectorField>> stagVectorStageSources_;

    // Private Member Functions

        // Return stage sources for given field type
        template<class FieldType>
        inline FastPtrList<FastPtrList<FieldType>>& stageSources();

        // Cycle stage sources for given field type
        template<class FieldType>
        inline void cycle();

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
    static autoPtr<RungeKuttaScheme> New(const fvMesh& fvMsh);

    // Construct from mesh and sizes
    RungeKuttaScheme
    (
        const fvMesh& fvMsh,
        const label nStages,
        const label nSteps = 1
    );

    // Destructor
    virtual ~RungeKuttaScheme();

    // Member Functions

        // Generate and return a list of fields to store the intermediate stage
        // sources or solutions, with given field name. The list stores all but
        // the last stage in each step.
        template<class Type, class MeshType>
        FastPtrList<meshField<Type,MeshType>>& newStageList(const word name);

        // Helper function for field types
        template<class FieldType>
        inline FastPtrList<FieldType>& newStageList(const word name)
        {
            return
                newStageList
                <
                    typename FieldType::dataType,
                    typename FieldType::meshType
                >(name);
        }

        // Access

            // Return the mesh
            inline const fvMesh& fvMsh() const { return fvMsh_; }

            // Return the Runge-Kutta dictionary
            inline const dictionary& dict() const { return dict_; }

            // Current stage level
            inline label stage() const { return stage_; }

            // Return the number of stages
            inline label nStages() const { return nStages_; }

            // Number of time steps
            inline label nSteps() const { return nSteps_; }

            // Current diagonal A coefficient
            inline scalar A() const;

            // Current diagonal B coefficient
            inline scalar B() const;

            // Current C coefficient
            inline scalar C() const;

            // Return the Butcher tableau for source A
            virtual const scalarRectangularMatrix& a() const { return a_; }

            // Return the Butcher tableau for source B
            virtual const scalarRectangularMatrix& b() const { return b_; }

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

        // Solution control

            // Loop over stages
            inline Switch loop();

            // Sum the given list of fields for the current stage, weighted by
            // the A Butcher tableau
            template<class Type, class MeshType>
            tmp<meshField<Type,MeshType>> stageSumA
            (
                const FastPtrList<meshField<Type,MeshType>>& list
            ) const;

            // Sum the given list of fields for the current stage, weighted by
            // the B Butcher tableau
            template<class Type, class MeshType>
            tmp<meshField<Type,MeshType>> stageSumB
            (
                const FastPtrList<meshField<Type,MeshType>>& list
            ) const;

            // Sum the two given lists of fields for the current stage, weighted
            // by the A and B Butcher tableaus, respectively
            template<class Type, class MeshType>
            tmp<meshField<Type,MeshType>> stageSum
            (
                const FastPtrList<meshField<Type,MeshType>>& listA,
                const FastPtrList<meshField<Type,MeshType>>& listB
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
