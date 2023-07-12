/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
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
    Foam::
    blendingLayerVelFvPatchVectorField

Description
    Reads lists of boundary values with spatial and temporal variation based 
    on mesoscale simulations

    Written by Aytac Kubilay    

    Example usage:
        myPatchName
        {
            type            blendingLayerVel;
            inputTimeStep   3600;
            value           uniform (1 1 0);
        }

SourceFiles
    blendingLayerVelFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef blendingLayerVelFvPatchVectorField_H
#define blendingLayerVelFvPatchVectorField_H

#include "fvPatchFields.H"
#include "fixedValueFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
    Class blendingLayerVelFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class blendingLayerVelFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{        
    // Private data
        
        //mesoscale timestep interval
        label inputTimeStep;

        //mesoscale target field
        vectorField Utarget_Field;                     

public:

    //- Runtime type information
    TypeName("blendingLayerVel");


    // Constructors

        //- Construct from patch and internal field
        blendingLayerVelFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        blendingLayerVelFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  turbulentTemperatureCoupledBaffleMixedFvPatchScalarField onto a
        //  new patch
        blendingLayerVelFvPatchVectorField
        (
            const
            blendingLayerVelFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );
        
        //- Construct as copy
        blendingLayerVelFvPatchVectorField
        (
            const blendingLayerVelFvPatchVectorField&
        );        

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new blendingLayerVelFvPatchVectorField
                (
                    *this
                )
            );
        }

        //- Construct as copy setting internal field reference
        blendingLayerVelFvPatchVectorField
        (
            const blendingLayerVelFvPatchVectorField&,
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
                new blendingLayerVelFvPatchVectorField
                (
                    *this,
                    iF
                )
            );
        }


    // Member functions

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
