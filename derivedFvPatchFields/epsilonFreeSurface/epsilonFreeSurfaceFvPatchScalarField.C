/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2018 OpenFOAM Foundation
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

#include "epsilonFreeSurfaceFvPatchScalarField.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "wallFvPatch.H"
#include "kEpsilon.H" 
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void epsilonFreeSurfaceFvPatchScalarField::writeLocalEntries(Ostream& os) const
{
    os.writeKeyword("up") << up_ << token::END_STATEMENT << nl;
}



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

epsilonFreeSurfaceFvPatchScalarField::epsilonFreeSurfaceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(p, iF),
    up_(0, 0, 1)
{
}


epsilonFreeSurfaceFvPatchScalarField::epsilonFreeSurfaceFvPatchScalarField
(
    const epsilonFreeSurfaceFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<scalar>(ptf, p, iF, mapper),
    up_(ptf.up_)
{
}


epsilonFreeSurfaceFvPatchScalarField::epsilonFreeSurfaceFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<scalar>(p, iF, dict),
	up_(dict.lookup("up"))
{
}


epsilonFreeSurfaceFvPatchScalarField::epsilonFreeSurfaceFvPatchScalarField
(
    const epsilonFreeSurfaceFvPatchScalarField& ptf
)
:
    fixedValueFvPatchField<scalar>(ptf),
    up_(ptf.up_)
{
}


epsilonFreeSurfaceFvPatchScalarField::epsilonFreeSurfaceFvPatchScalarField
(
    const epsilonFreeSurfaceFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchField<scalar>(ptf, iF),
    up_(ptf.up_)
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void epsilonFreeSurfaceFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

	// Get range and orientation
	boundBox bb(patch().patch().localPoints(), true);

	// Height of free surface
	scalar height = bb.max() & up_;

	// From v2WallFunctionfvPatchScalarField.C
    const turbulenceModel& turbModel = db().lookupObject<turbulenceModel>
    (
        IOobject::groupName
        (
            turbulenceModel::propertiesName,
            internalField().group()
        )
    );

	const tmp<volScalarField> tk = turbModel.k();
    const volScalarField& k = tk();


    scalarField& epsilon = *this;

    // Set epsilon patch values
    forAll(epsilon, facei)
    {
		label celli = patch().faceCells()[facei];
		epsilon[facei] = pow(k[celli], 1.5)/(0.43*height);
	}

	fixedValueFvPatchField<scalar>::updateCoeffs();

    // TODO: perform averaging for cells sharing more than one boundary face
}


void epsilonFreeSurfaceFvPatchScalarField::evaluate
(
    const Pstream::commsTypes commsType
)
{
    fixedValueFvPatchField<scalar>::evaluate(commsType);
}


void epsilonFreeSurfaceFvPatchScalarField::write(Ostream& os) const
{
    writeLocalEntries(os);
    fixedValueFvPatchField<scalar>::write(os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    epsilonFreeSurfaceFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
