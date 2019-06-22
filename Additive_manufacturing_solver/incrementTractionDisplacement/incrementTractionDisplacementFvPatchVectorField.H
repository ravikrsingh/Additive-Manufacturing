/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::tractionDisplacementFvPatchVectorField

Description
    Fixed traction boundary condition for the standard linear elastic, fixed
    coefficient displacement equation.

SourceFiles
    tractionDisplacementFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef incrementTractionDisplacementFvPatchVectorField_H
#define incrementTractionDisplacementFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class incrementTractionDisplacementFvPatch Declaration
\*---------------------------------------------------------------------------*/

class incrementTractionDisplacementFvPatchVectorField
:
    public fixedGradientFvPatchVectorField
{

    // Private Data

        vectorField incrementTraction_;
        scalarField incrementPressure_;


public:

    //- Runtime type information
    TypeName("incrementTractionDisplacement");


    // Constructors

        //- Construct from patch and internal field
        incrementTractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        incrementTractionDisplacementFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  incrementTractionDisplacementFvPatchVectorField onto a new patch
        incrementTractionDisplacementFvPatchVectorField
        (
            const incrementTractionDisplacementFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        incrementTractionDisplacementFvPatchVectorField
        (
            const incrementTractionDisplacementFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new incrementTractionDisplacementFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        incrementTractionDisplacementFvPatchVectorField
        (
            const incrementTractionDisplacementFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new incrementTractionDisplacementFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

        // Access

            virtual const vectorField& incrementTraction() const
            {
                return incrementTraction_;
            }

            virtual vectorField& incrementTraction()
            {
                return incrementTraction_;
            }

            virtual const scalarField& incrementPressure() const
            {
                return incrementPressure_;
            }

            virtual  scalarField& incrementPressure()
            {
                return incrementPressure_;
            }


        // Mapping functions

            //- Map (and resize as needed) from self given a mapping object
            virtual void autoMap
            (
                const fvPatchFieldMapper&
            );

            //- Reverse map the given fvPatchField onto this fvPatchField
            virtual void rmap
            (
                const fvPatchVectorField&,
                const labelList&
            );


        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
