/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2014 OpenFOAM Foundation
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

#include "simpleGrass.H"
#include "fvmLaplacian.H"
#include "fvmSup.H"
#include "constants.H"
#include "addToRunTimeSelectionTable.H"

#include "mixedFvPatchFields.H"
#include "mappedPatchBase.H"

using namespace Foam::constant;

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    namespace grass
    {
        defineTypeNameAndDebug(simpleGrass, 0);
        addToGrassRunTimeSelectionTables(simpleGrass);
    }
}

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

void Foam::grass::simpleGrass::initialise()
{
    List<word> grassPatches(coeffs_.lookup("grassPatches"));   
    LAI_ = coeffs_.lookupOrDefault("LAI", 2);
    nEvapSides_ = coeffs_.lookupOrDefault("nEvapSides", 1);
    beta_ = coeffs_.lookupOrDefault("beta", 0.78);

    label count = 0;
    forAll(grassPatches, i)
    {
        grassPatchID = mesh_.boundaryMesh().findPatchID(grassPatches[i]);
        if (grassPatchID < 0)
        {
            FatalErrorInFunction
                << "Grass patch named " << grassPatches[i] << " not found." << nl
                << abort(FatalError);
        }
        else
        {
            selectedPatches_[count] = grassPatchID;
            count++;
        }
    }
    selectedPatches_.resize(count--);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::grass::simpleGrass::simpleGrass(const volScalarField& T)
:
    grassModel(typeName, T),
    Tg_
    (
        IOobject
        (
            "Tg",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Tg", dimTemperature, 300),
		T.boundaryField().types()
    ),
    T
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    w
    (
        IOobject
        (
            "w",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    qs
    (
        IOobject
        (
            "qs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    qr
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    selectedPatches_(mesh_.boundary().size(), -1)
{
    initialise();
}


Foam::grass::simpleGrass::simpleGrass(const dictionary& dict, const volScalarField& T)
:
    grassModel(typeName, dict, T),
    Tg_
    (
        IOobject
        (
            "Tg",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("Tg", dimTemperature, 300),
		T.boundaryField().types()
    ),
    T
    (
        IOobject
        (
            "T",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    w
    (
        IOobject
        (
            "w",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    qs
    (
        IOobject
        (
            "qs",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    qr
    (
        IOobject
        (
            "qr",
            mesh_.time().timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh_
    ),
    selectedPatches_(mesh_.boundary().size(), -1)
{
    initialise();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::grass::simpleGrass::~simpleGrass()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::grass::simpleGrass::read()
{
    if (grassModel::read())
    {
        // nothing to read

        return true;
    }
    else
    {
        return false;
    }
}


void Foam::grass::simpleGrass::calculate()
{	

	forAll(selectedPatches_, i)
    {
		label patchID = selectedPatches_[i];
		const fvPatch& thisPatch = mesh_.boundary()[patchID];

        volScalarField::Boundary& TgBf = Tg_.boundaryFieldRef();
        scalarField& TgPatch = TgBf[patchID];
        volScalarField::Boundary& qsBf = qs.boundaryFieldRef();
        scalarField& qsPatch = qsBf[patchID];
        volScalarField::Boundary& qrBf = qr.boundaryFieldRef();
        scalarField& qrPatch = qrBf[patchID];

        scalarField TPatchInternal = thisPatch.patchInternalField(T);
        scalarField wPatchInternal = thisPatch.patchInternalField(w);
        scalarField TPatch = thisPatch.lookupPatchField<volScalarField, scalar>("T");

		calc_leafTemp(TgPatch, TPatch, qsPatch, qrPatch, TPatchInternal, wPatchInternal);
	}
}

void Foam::grass::simpleGrass::calc_leafTemp
(
	scalarField& leafTemp,
	const scalarField& Ts,  
	scalarField& Qs, 
	scalarField& Qr,
	const scalarField& Tc,
	const scalarField& wc
)
{
	const scalar p_ = 101325; // air pressure - Pa
	const scalar rhoa = 1.225; // density of air - kg/m3
	const scalar cpa = 1003.5; // specific heat of air at constant pressure - J/(kgK)
	scalar rs = 200; //stomatal resistance - s/m
	scalar ra = 100; //aerodynamic resistance - s/m

    scalarField Qs_abs = Qs - Qs*exp(beta_*LAI_);
    Qs = Qs*exp(beta_*LAI_);
    scalarField Qr_abs = Qr;
    Qr = 0;

	label maxIter = 500;
	for (label i=1; i<=maxIter; i++)
	{
		scalarField evsat = exp( - 5.8002206e3/leafTemp // saturated vapor pressure pws - ASHRAE 1.2
	                + 1.3914993
	                - 4.8640239e-2*leafTemp
	                + 4.1764768e-5*pow(leafTemp,2)
	                - 1.4452093e-8*pow(leafTemp,3)
	                + 6.5459673*log(leafTemp) );

		scalarField rhosat = evsat / (461.5*leafTemp); // saturated density of water vapour

		scalarField wsat = 0.621945*(evsat/(p_-evsat)); // saturated specific humidity, ASHRAE 1, eq.23

		///////////calculate transpiration rate
		scalarField E = pos(Qs-SMALL)*nEvapSides_*LAI_*rhoa*(wsat-wc)/(rs+ra);
		//no transpiration at night when Qs is not >0
		//////////////////////////////////////

		scalar lambda = 2500000; // latent heat of vaporization of water J/Kg
		scalarField Qlat = lambda*E; //latent heat flux

        scalarField Qr_Surface = 6*(Ts-leafTemp); //thermal radiation between grass and surface - Malys et al 2014
		scalarField leafTemp_new = Tc + (Qr_abs + Qr_Surface + Qs_abs - Qlat)*(ra/(rhoa*cpa*2*LAI_));

        // info
        Info << " max leaf temp tl=" << gMax(leafTemp_new)
             << " K, iteration i="   << i << endl;

        // Check rel. L-infinity error
        scalar maxError = gMax(mag(leafTemp_new-leafTemp));
        scalar maxRelError = maxError/gMax(mag(leafTemp_new));  

        // update leaf temp.
        //Tl_.internalField() = new_Tl.internalField();
        //Tl_.boundaryField() = new_Tl.boundaryField();
        leafTemp = 0.5*leafTemp+0.5*leafTemp_new;  

        // convergence check
        if (maxRelError < 1e-8)
             break;                 
	}
}


// ************************************************************************* //
