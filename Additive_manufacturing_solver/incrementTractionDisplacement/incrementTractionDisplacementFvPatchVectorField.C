/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "incrementTractionDisplacementFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "fvCFD.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

incrementTractionDisplacementFvPatchVectorField::
incrementTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(p, iF),
    incrementTraction_(p.size(), Zero),
    incrementPressure_(p.size(), 0.0)
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = Zero;
}


incrementTractionDisplacementFvPatchVectorField::
incrementTractionDisplacementFvPatchVectorField
(
    const incrementTractionDisplacementFvPatchVectorField& tdpvf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedGradientFvPatchVectorField(tdpvf, p, iF, mapper),
    incrementTraction_(tdpvf.incrementTraction_, mapper),
    incrementPressure_(tdpvf.incrementPressure_, mapper)
{}


incrementTractionDisplacementFvPatchVectorField::
incrementTractionDisplacementFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchVectorField(p, iF),
    incrementTraction_("incrementTraction", dict, p.size()),
    incrementPressure_("incrementPressure", dict, p.size())
{
    fvPatchVectorField::operator=(patchInternalField());
    gradient() = Zero;
}


incrementTractionDisplacementFvPatchVectorField::
incrementTractionDisplacementFvPatchVectorField
(
    const incrementTractionDisplacementFvPatchVectorField& tdpvf
)
:
    fixedGradientFvPatchVectorField(tdpvf),
    incrementTraction_(tdpvf.incrementTraction_),
    incrementPressure_(tdpvf.incrementPressure_)
{}


incrementTractionDisplacementFvPatchVectorField::
incrementTractionDisplacementFvPatchVectorField
(
    const incrementTractionDisplacementFvPatchVectorField& tdpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedGradientFvPatchVectorField(tdpvf, iF),
    incrementTraction_(tdpvf.incrementTraction_),
    incrementPressure_(tdpvf.incrementPressure_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void incrementTractionDisplacementFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedGradientFvPatchVectorField::autoMap(m);
    incrementTraction_.autoMap(m);
    incrementPressure_.autoMap(m);
}


void incrementTractionDisplacementFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedGradientFvPatchVectorField::rmap(ptf, addr);

    const incrementTractionDisplacementFvPatchVectorField& dmptf =
        refCast<const incrementTractionDisplacementFvPatchVectorField>(ptf);

    incrementTraction_.rmap(dmptf.incrementTraction_, addr);
    incrementPressure_.rmap(dmptf.incrementPressure_, addr);
}


void incrementTractionDisplacementFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

//     const dictionary& transportProperties =
//         db().lookupObject<IOdictionary>("transportProperties");


    vectorField n(patch().nf());

        
            const fvPatchField<scalar>&  muE=
            patch().lookupPatchField<volScalarField, scalar>("muE");
              
            const fvPatchField<scalar>&  lambdaE=
            patch().lookupPatchField<volScalarField, scalar>("lambdaE");
//     const fvPatchField<symmTensor>& sigmaD =
//         patch().lookupPatchField<volSymmTensorField, symmTensor>("sigmaD");
    const fvPatchField<tensor>& graddD =
    patch().lookupPatchField<volTensorField, vector>("graddD");
    const fvPatchField<symmTensor>& deps_p =
    patch().lookupPatchField<volSymmTensorField, tensor>("deps_p");
    const fvPatchField<symmTensor>& deps_eigen =
    patch().lookupPatchField<volSymmTensorField, symmTensor>("deps_eigen");
    gradient() =
    (
        (incrementTraction_ - incrementPressure_*n)
        - (n & (muE*graddD.T() - (muE + lambdaE)*graddD))
        - n*tr(graddD)*lambdaE
        + (n & (2*muE*deps_p))
        + n*tr(deps_p)*lambdaE
        + (n & (2*muE*deps_eigen))
        + n*tr(deps_eigen)*lambdaE
    )/(2.0*muE + lambdaE);
    

    fixedGradientFvPatchVectorField::updateCoeffs();
}


void incrementTractionDisplacementFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    incrementTraction_.writeEntry("incrementTraction", os);
    incrementPressure_.writeEntry("incrementPressure", os);
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    incrementTractionDisplacementFvPatchVectorField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
