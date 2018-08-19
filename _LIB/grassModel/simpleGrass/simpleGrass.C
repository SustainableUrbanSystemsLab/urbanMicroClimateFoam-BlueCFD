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
    nEvapSides_ = coeffs_.lookupOrDefault("nEvapSides", 1);
    Cd_ = coeffs_.lookupOrDefault("Cd", 0.2);
    beta_ = coeffs_.lookupOrDefault("beta", 0.78);
    LAI_ = coeffs_.lookupOrDefault("LAI", 2);
    LAD_ = coeffs_.lookupOrDefault("LAD", 20);
    p_ = coeffs_.lookupOrDefault("p", 101325);
    rhoa = coeffs_.lookupOrDefault("rhoa", 1.225);
    cpa = coeffs_.lookupOrDefault("cpa", 1003.5);
    rs = coeffs_.lookupOrDefault("rs", 200);
    ra = coeffs_.lookupOrDefault("ra", 100);

    List<word> grassPatches(coeffs_.lookup("grassPatches"));  

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
    Tl_
    (
        IOobject
        (
            "Tl",
            mesh_.time().timeName(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        -pos(T)
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


void Foam::grass::simpleGrass::calculate
(
    const volScalarField& T_, 
    const volScalarField& w_,
    volScalarField& Sh_,
    volScalarField& Sw_
)
{	

	forAll(selectedPatches_, patchi)
    {
		label patchID = selectedPatches_[patchi];
		const fvPatch& thisPatch = mesh_.boundary()[patchID];

        const labelUList& patchInternalLabels = thisPatch.faceCells();
        scalarField leafTemp = thisPatch.patchInternalField(Tl_);
        scalarField TPatchInternal = thisPatch.patchInternalField(T_);
        scalarField wPatchInternal = thisPatch.patchInternalField(w_);

        scalarField TPatch = thisPatch.lookupPatchField<volScalarField, scalar>("T");
        scalarField qsPatch = thisPatch.lookupPatchField<volScalarField, scalar>("qs");
        scalarField qrPatch = thisPatch.lookupPatchField<volScalarField, scalar>("qr");

        scalarField transpiration(scalarField(mesh_.nCells(),0.0));

        if (gMin(leafTemp) < 0)
        {
            leafTemp = TPatchInternal; //initialize if necessary
        }      

		calc_leafTemp(leafTemp, transpiration, TPatch, qsPatch, qrPatch, TPatchInternal, wPatchInternal);

        scalarField& Tl_Internal = Tl_.ref();
        volScalarField::Boundary& Tl_Bf = Tl_.boundaryFieldRef();
        scalarField& Tl_p = Tl_Bf[patchID];

        forAll(patchInternalLabels, i)
        {
            Tl_Internal[patchInternalLabels[i]] = leafTemp[i];
            Tl_p[i] = leafTemp[i];
            Sh_[patchInternalLabels[i]] = 2.0*rhoa*cpa*LAD_*(leafTemp[i]-TPatchInternal[i])/ra;
            Sw_[patchInternalLabels[i]] = transpiration[i]*LAD_;
        }

	}
}

void Foam::grass::simpleGrass::calc_leafTemp
(
	scalarField& leafTemp,
	scalarField& E,
	const scalarField& Ts,  
	const scalarField& Qs, 
	const scalarField& Qr,
	const scalarField& Tc,
	const scalarField& wc
)
{

    scalarField Qs_abs = Qs - Qs*exp(-beta_*LAI_);
//    Qs = Qs*exp(-beta_*LAI_);
    scalarField Qr_abs = Qr;
//    Qr = 0;

	label maxIter = 500;
	for (label i=1; i<=maxIter; i++)
	{
		scalarField evsat = exp( - 5.8002206e3/leafTemp // saturated vapor pressure pws - ASHRAE 1.2
	                + 1.3914993
	                - 4.8640239e-2*leafTemp
	                + 4.1764768e-5*pow(leafTemp,2)
	                - 1.4452093e-8*pow(leafTemp,3)
	                + 6.5459673*log(leafTemp) );

		scalarField wsat = 0.621945*(evsat/(p_-evsat)); // saturated specific humidity, ASHRAE 1, eq.23

		///////////calculate transpiration rate
		E = pos(Qs-SMALL)*nEvapSides_*rhoa*(wsat-wc)/(rs+ra);
		//no transpiration at night when Qs is not >0
		//////////////////////////////////////

		scalar lambda = 2500000; // latent heat of vaporization of water J/kg
		scalarField Qlat = lambda*E*LAI_; //latent heat flux

        scalarField Qr_Surface = 6*(Ts-leafTemp); //thermal radiation between grass and surface - Malys et al 2014
		scalarField leafTemp_new = Tc + (Qr_abs + Qr_Surface + Qs_abs - Qlat)*(ra/(rhoa*cpa*2*LAI_));

        // info
        Info << " max leaf temp Tl=" << gMax(leafTemp_new)
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

Foam::tmp<Foam::volScalarField> Foam::grass::simpleGrass::Cf() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Cf",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false
            ),
            dimensionedScalar("Cf", dimless/dimLength, Cd_*LAD_)*pos(Tl_)
        )
    );
}

// ************************************************************************* //
