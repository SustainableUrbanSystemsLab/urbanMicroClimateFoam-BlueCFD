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
    betaLW_ = coeffs_.lookupOrDefault("betaLW", 0.83);
    LAI_ = coeffs_.lookupOrDefault("LAI", 2.0);
    LAD_ = coeffs_.lookupOrDefault("LAD", 20.0);
    albedoSoil_ = coeffs_.lookupOrDefault("albedoSoil", 0.2366);
    emissivitySoil_ = coeffs_.lookupOrDefault("emissivitySoil", 0.95);

    p_ = 101325;
    rhoa = 1.225;
    cpa = 1003.5;
    Ra = 287.042;
    Rv = 461.524;

    rs = coeffs_.lookupOrDefault("rs", 200);
//    ra = coeffs_.lookupOrDefault("ra", 100);
    debug_ = coeffs_.lookupOrDefault("debug", 0);

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

Foam::scalarField Foam::grass::simpleGrass::calc_pvsat(const scalarField& T_)
{
     scalarField pvsat_ = exp( - 5.8002206e3/T_ // saturated vapor pressure pws - ASHRAE 1.2
            + 1.3914993
            - 4.8640239e-2*T_
            + 4.1764768e-5*pow(T_,2)
            - 1.4452093e-8*pow(T_,3)
            + 6.5459673*log(T_) );
     return pvsat_;
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
    const volVectorField& U_,
    volScalarField& Sh_,
    volScalarField& Sw_
)
{    

    forAll(selectedPatches_, patchi)
    {
        label patchID = selectedPatches_[patchi];
        const fvPatch& thisPatch = mesh_.boundary()[patchID];

        const labelUList& patchInternalLabels = thisPatch.faceCells();
        scalarField Tl = thisPatch.patchInternalField(Tl_);
        scalarField Tc = thisPatch.patchInternalField(T_);
        scalarField wc = thisPatch.patchInternalField(w_);
        vectorField Uc = thisPatch.patchInternalField(U_);

        // Get the coupling information from the mappedPatchBase
        const mappedPatchBase& mpp =
            refCast<const mappedPatchBase>(thisPatch.patch());
        const polyMesh& nbrMesh = mpp.sampleMesh();
        const label samplePatchI = mpp.samplePolyPatch().index();
        const fvPatch& nbrPatch =
            refCast<const fvMesh>(nbrMesh).boundary()[samplePatchI];
        scalarField Ts = nbrPatch.lookupPatchField<volScalarField, scalar>("Ts");
            mpp.distribute(Ts);
// read Ts from solid domain instead
//        scalarField Ts = thisPatch.lookupPatchField<volScalarField, scalar>("T");

        scalarField qs = thisPatch.lookupPatchField<volScalarField, scalar>("qs");
        scalarField qr = thisPatch.lookupPatchField<volScalarField, scalar>("qr");

        scalar C_ = 131.035; // proportionality factor
        scalar l_ = 0.1; // characteristic length of leaf
        scalarField magU = mag(Uc);
        magU = max(magU, SMALL);
        ra = C_*pow(l_/magU, 0.5); //aerodynamic resistance

        scalarField h_ch = (2*rhoa*cpa)/ra; //convective heat transfer coefficient
        scalarField h_cm = (rhoa*Ra)/(p_*Rv*(ra+rs)); //convective mass transfer coefficient

        //scalarField wsat = 0.621945*(pvsat/(p_-pvsat)); // saturated specific humidity, ASHRAE 1, eq.23
        scalarField pv = p_ * wc / (Ra/Rv + wc); // vapour pressure

        if (gMin(Tl) < 0)
        {
            Tl = Tc; //initialize if necessary
        }

        scalarField Qs_abs = qs*(1-exp(-beta_*LAI_))*(1+exp(-beta_*LAI_)*albedoSoil_);

        scalarField E(scalarField(Qs_abs.size(),0.0));

        ////calculate leaf temperature///////////
        label maxIter = 500;
        for (label i=1; i<=maxIter; i++)
        {
            scalarField pvsat = calc_pvsat(Tl); //saturation vapour pressure
            E = pos(Qs_abs-SMALL)*nEvapSides_*h_cm*(pvsat-pv); //initialize transpiration rate [kg/(m2s)]
            //scalarField E = pos(Qs-SMALL)*nEvapSides_*rhoa*(wsat-wc)/(rs+ra);
            //no transpiration at night when Qs_abs is not >0

            scalar lambda = 2500000; // latent heat of vaporization of water J/kg
            scalarField Qlat = lambda*E*LAI_; //latent heat flux

            scalarField Qr2surrounding = qr;
            scalarField Qr2substrate = 6*(Ts-Tl); //thermal radiation between grass and surface - Malys et al 2014

            scalarField Tl_new = Tc + (Qr2surrounding + Qr2substrate + Qs_abs - Qlat)/ (h_ch*LAI_);

            // info
            Info << " max leaf temp Tl=" << gMax(Tl_new)
                 << " K, iteration i="   << i << endl;

            // Check rel. L-infinity error
            scalar maxError = gMax(mag(Tl_new-Tl));
            scalar maxRelError = maxError/gMax(mag(Tl_new));  

            // convergence check
            if ((maxRelError < 1e-8) && (maxError < 1e-8))
            {
                if(debug_)
                {
                    scalarField Qsen = h_ch*(Tc-Tl)*LAI_;
                    Info << " Qs_abs: " << gSum(thisPatch.magSf()*Qs_abs)/gSum(thisPatch.magSf()) << endl;
                    Info << " Qlat: " << gSum(thisPatch.magSf()*-Qlat)/gSum(thisPatch.magSf()) << endl;
                    Info << " Qsen: " << gSum(thisPatch.magSf()*Qsen)/gSum(thisPatch.magSf()) << endl;
                    Info << " Qr2surrounding: " << gSum(thisPatch.magSf()*Qr2surrounding)/gSum(thisPatch.magSf()) << endl;
                    Info << " Qr2substrate: " << gSum(thisPatch.magSf()*Qr2substrate)/gSum(thisPatch.magSf()) << endl;
                }
                break; 
            }
            else
            {
                Tl = 0.5*Tl+0.5*Tl_new; // update leaf temp. 
            }                    
        }
        /////////////////////////////////////////

        ////update fields////////////////////////
        scalarField& Tl_Internal = Tl_.ref();
        volScalarField::Boundary& Tl_Bf = Tl_.boundaryFieldRef();
        scalarField& Tl_p = Tl_Bf[patchID];

        forAll(patchInternalLabels, i)
        {
            Tl_Internal[patchInternalLabels[i]] = Tl[i];
            Tl_p[i] = Tl[i];
            Sh_[patchInternalLabels[i]] = LAD_*h_ch[i]*(Tl[i]-Tc[i]);
            Sw_[patchInternalLabels[i]] = LAD_*E[i];
        }
        /////////////////////////////////////////
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

