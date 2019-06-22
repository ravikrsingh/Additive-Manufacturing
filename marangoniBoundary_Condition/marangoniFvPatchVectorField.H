/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2007 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Class
    marangoniFvPatchVectorField

Description

SourceFiles
    marangoniFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef marangoniFvPatchVectorField_H
#define marangoniFvPatchVectorField_H
#include "fvPatchFields.H"
#include "transformFvPatchFields.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                      Class marangoniFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class marangoniFvPatchVectorField
:
    public transformFvPatchField<vector>
{
    // Private data

    // proportionality factor btween thermal and momntum gradient
    scalar marangonicoeff_;

public:

    //- Runtime type information
    TypeName("marangoni");

    // Constructors

        //- Construct from patch and internal field
        marangoniFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        marangoniFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given marangoniFvPatchVectorField
        //  onto a new patch
        marangoniFvPatchVectorField
        (
            const marangoniFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        marangoniFvPatchVectorField
        (
            const marangoniFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField > clone() const
        {
            return tmp<fvPatchVectorField >
            (
                new marangoniFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        marangoniFvPatchVectorField
        (
            const marangoniFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField > clone(const DimensionedField<vector, volMesh>& iF) const
        {
            return tmp<fvPatchVectorField>
            (
                new marangoniFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

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


        // Return defining fields

    // virtual vectorField& fixedValue()
    //       {
    //          return fixedValue_;
    //      }

    //      virtual const vectorField& fixedValue() const
    //      {
    //          return fixedValue_;
    //      }


        // Evaluation functions

            //- Return gradient at boundary
            virtual tmp<vectorField > snGrad() const;

            //- Evaluate the patch field
            virtual void evaluate();

            //- Return face-gradient transform diagonal
            virtual tmp<vectorField > snGradTransformDiag() const;


        //- Write
        virtual void write(Ostream&) const;


    // Member operators

    virtual void operator=(const UList<vector>&) {}

    virtual void operator=(const fvPatchVectorField&) {}
      virtual void operator+=(const fvPatchVectorField&) {}
      virtual void operator-=(const fvPatchVectorField&) {}
      virtual void operator*=(const fvPatchField<scalar>&) {}
      virtual void operator/=(const fvPatchField<scalar>&) {}
      virtual void operator+=(const vectorField&) {}
      virtual void operator-=(const vectorField&) {}

      virtual void operator*=(const Field<scalar>&) {}
      virtual void operator/=(const Field<scalar>&) {}


            virtual void operator=(const vector&) {}
      virtual void operator+=(const vector&) {}
      virtual void operator-=(const vector&) {}
      virtual void operator*=(const scalar) {}
      virtual void operator/=(const scalar) {}
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

//#ifdef NoRepository
//#   include "marangoniFvPatchVectorField.C"
//#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

