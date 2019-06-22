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

\*---------------------------------------------------------------------------*/

#include "marangoniFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchField.H"
#include "volFields.H"
#include "transformFvPatchFields.H"
#include "symmTransformField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

marangoniFvPatchVectorField::marangoniFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchVectorField(p, iF),
    marangonicoeff_(0)
{
    //this->checkVolField();
}


marangoniFvPatchVectorField::marangoniFvPatchVectorField
(
    const marangoniFvPatchVectorField& ptf,
    const fvPatch& p,
   const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    transformFvPatchVectorField(ptf, p, iF, mapper),
    marangonicoeff_(ptf.marangonicoeff_)
{
    //this->checkVolField();
}


marangoniFvPatchVectorField::marangoniFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    transformFvPatchVectorField(p, iF)
{
    //this->checkVolField();
    marangonicoeff_=readScalar(dict.lookup("marangonicoeff"));
    evaluate();
}


marangoniFvPatchVectorField::marangoniFvPatchVectorField
(
    const marangoniFvPatchVectorField& ptf
)
:
    transformFvPatchVectorField(ptf),
    marangonicoeff_(ptf.marangonicoeff_)
{
    //this->checkVolField();
}


marangoniFvPatchVectorField::marangoniFvPatchVectorField
(
    const marangoniFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    transformFvPatchVectorField(ptf, iF),
    marangonicoeff_(ptf. marangonicoeff_)
{
    //this->checkVolField();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Map from self
void marangoniFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    vectorField::autoMap(m);
}


// Reverse-map the given fvPatchField onto this fvPatchField
void marangoniFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    transformFvPatchVectorField::rmap(ptf, addr);

    const marangoniFvPatchVectorField& dmptf =
        refCast<const marangoniFvPatchVectorField >(ptf);

    marangonicoeff_=dmptf.marangonicoeff_;
}


// Return gradient at boundary
tmp<vectorField > marangoniFvPatchVectorField::snGrad() const
{
    //Info << "entering  marangoniFvPatchVectorField::snGrad()" << endl;
    vectorField nHat = this->patch().nf();
    vectorField pif = this->patchInternalField();
    vectorField result;
    if(!db().foundObject<vectorField>("gradT")) {
 	Info << " marangoniFvPatchVectorField::snGrad(): object gradT not found!" << endl;
 	Foam::error theError(" marangoniFvPatchVectorField::evaluate(): object gradT not found!");
	return transform(I - sqr(nHat),pif);
    }
    else {
	
 	fvPatchField<vector> gradT =this->patch().lookupPatchField<volVectorField, vector>("gradT");
 	vectorField  gradT_internal = gradT.patchInternalField();
 	vectorField gradTplane= transform(I - sqr(nHat),gradT_internal);
 	vectorField pifplane= transform(I - sqr(nHat),pif);
	result=pifplane+marangonicoeff_*gradTplane/this->patch().deltaCoeffs();
	return (result-pif)*this->patch().deltaCoeffs();
    }
    //Info << "leaving  marangoniFvPatchVectorField::snGrad()" << endl;    
    return result;
}


// Evaluate the field on the patch
void marangoniFvPatchVectorField::evaluate()
{
    //Info << "entering  marangoniFvPatchVectorField::evaluate()" << endl;
    if (!this->updated())
    {
	//Info << "marangoniFvPatchVectorField::evaluate(): calling updatecoeffs" << endl;
        this->updateCoeffs();
    }

     vectorField nHat = this->patch().nf();
     vectorField pif = this->patchInternalField();
     scalarField deltas=this->patch().deltaCoeffs();

     if(!db().foundObject<vectorField>("gradT")) {
 	Info << " marangoniFvPatchVectorField::snGrad(): object gradT not found!" << endl;
 	Foam::error theError(" marangoniFvPatchVectorField::evaluate(): object gradT not found!");
 	vectorField::operator=
 	    (
 	      	transform(I - sqr(nHat),pif)
 	     );
 	// theError.exit();
     }
     else {

 	fvPatchField<vector> gradT =this->patch().lookupPatchField<volVectorField, vector>("gradT");
 	vectorField  gradT_internal = gradT.patchInternalField();
 	vectorField gradTplane= transform(I - sqr(nHat),gradT_internal);
 	vectorField pifplane= transform(I - sqr(nHat),pif);
	vectorField result=pifplane+marangonicoeff_*gradTplane/deltas;

 	vectorField::operator=
 	    (
	     result
 	     );
     }

    transformFvPatchVectorField::evaluate();
}


// Return defining fields
tmp<vectorField > marangoniFvPatchVectorField::snGradTransformDiag() const
{
    vectorField nHat = this->patch().nf();
    vectorField diag(nHat.size());

    diag.replace(vector::X, mag(nHat.component(vector::X)));
    diag.replace(vector::Y, mag(nHat.component(vector::Y)));
    diag.replace(vector::Z, mag(nHat.component(vector::Z)));

    return transformFieldMask<vector>(pow<vector, pTraits<vector>::rank>(diag));
}


// Write
void marangoniFvPatchVectorField::write(Ostream& os) const
{
    transformFvPatchVectorField::write(os);
    os.writeKeyword("marangonicoeff") << marangonicoeff_ << token::END_STATEMENT << nl;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField(fvPatchVectorField, marangoniFvPatchVectorField);


} // End namespace Foam

// ************************************************************************* //
