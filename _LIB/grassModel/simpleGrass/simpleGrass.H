/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::grass::simpleGrass

Description
    Works well for combustion applications where optical thickness, tau is
    large, i.e. tau = a*L > 3 (L = distance between objects)

    Assumes
     - all surfaces are diffuse
     - tends to over predict radiative fluxes from sources/sinks
       *** SOURCES NOT CURRENTLY INCLUDED ***

SourceFiles
    simpleGrass.C

\*---------------------------------------------------------------------------*/

#ifndef grassModelsimpleGrass_H
#define grassModelsimpleGrass_H

#include "grassModel.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace grass
{

/*---------------------------------------------------------------------------*\
                           Class simpleGrass Declaration
\*---------------------------------------------------------------------------*/

class simpleGrass
:
    public grassModel
{
    // Private data

        //- leaf temperature [K]
        volScalarField Tl_;

        //- Selected patches
        labelList selectedPatches_;

        scalar LAI_; // leaf area index
        scalar LAD_; // leaf area density
        dimensionedScalar Cf_; // leaf drag coefficient
        scalar nEvapSides_; // number of sides, leaf evaporates from
        scalar beta_; // extinction coefficient for short-wave radiation
		scalar p_; // air pressure - Pa
		scalar rhoa; // density of air - kg/m3
		scalar cpa; // specific heat of air at constant pressure - J/(kgK)
		scalar rs; //stomatal resistance - s/m
		scalar ra; //aerodynamic resistance - s/m

        //- Grass patch ID
        label grassPatchID;

    // Private Member Functions

        //- Initialise
        void initialise();

        //- Disallow default bitwise copy construct
        simpleGrass(const simpleGrass&);

        //- Disallow default bitwise assignment
        void operator=(const simpleGrass&);


public:

    //- Runtime type information
    TypeName("simpleGrass");


    // Constructors

        //- Construct from components
        simpleGrass(const volScalarField& T);


    //- Destructor
    virtual ~simpleGrass();


    // Member functions

        // Edit

            //- Solve grass equation(s)
            void calculate
            (
                const volScalarField& T_, 
                const volScalarField& w_,
                volScalarField& Sh_,
                volScalarField& Sw_
            );

            //- calculate leaf temperature
            void calc_leafTemp
			(
				scalarField& leafTemp, 
                scalarField& transpiration,
          	    const scalarField& Ts,  
				const scalarField& Qs, 
				const scalarField& Qr,
				const scalarField& Tc,
				const scalarField& wc
			);

            //- Read grass properties dictionary
            bool read();

            //- return leaf drag coefficient
            virtual tmp<volScalarField> Cf() const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace grass
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
