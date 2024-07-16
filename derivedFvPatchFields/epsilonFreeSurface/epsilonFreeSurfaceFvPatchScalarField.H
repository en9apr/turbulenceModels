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

Class
    Foam::RASModels::epsilonFreeSurfaceFvPatchScalarField

Description
    This boundary condition provides a turbulence damping function, f, wall
    function condition for low- and high Reynolds number, turbulent flow cases

    The model operates in two modes, based on the computed laminar-to-turbulent
    switch-over y+ value derived from kappa and E.

Usage
    \table
        Property     | Description             | Required    | Default value
        Cmu          | model coefficient       | no          | 0.09
        kappa        | Von Karman constant     | no          | 0.41
        E            | model coefficient       | no          | 9.8
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            epsilonFreeSurface;
    }
    \endverbatim

See also
    Foam::fixedValueFvPatchField

SourceFiles
    epsilonFreeSurfaceFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef epsilonFreeSurfaceFvPatchScalarField_H
#define epsilonFreeSurfaceFvPatchScalarField_H

#include "fixedValueFvPatchField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
               Class epsilonFreeSurfaceFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class epsilonFreeSurfaceFvPatchScalarField
:
    public fixedValueFvPatchField<scalar>
{
protected:

    // Protected data

		//- Direction of up
		vector up_;

    // Protected Member Functions

        //- Write local wall function variables
        virtual void writeLocalEntries(Ostream&) const;


public:

    //- Runtime type information
    TypeName("epsilonFreeSurface");


    // Constructors

        //- Construct from patch and internal field
        epsilonFreeSurfaceFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        epsilonFreeSurfaceFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given epsilonFreeSurfaceFvPatchScalarField
        //  onto a new patch
        epsilonFreeSurfaceFvPatchScalarField
        (
            const epsilonFreeSurfaceFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        epsilonFreeSurfaceFvPatchScalarField
        (
            const epsilonFreeSurfaceFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new epsilonFreeSurfaceFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        epsilonFreeSurfaceFvPatchScalarField
        (
            const epsilonFreeSurfaceFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new epsilonFreeSurfaceFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        // Mapping functions

			//- Return up direction
			vector& up()
			{
				return up_;
			}

        // Evaluation functions

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Evaluate the patchField
            virtual void evaluate(const Pstream::commsTypes);


        // I-O

            //- Write
            virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
