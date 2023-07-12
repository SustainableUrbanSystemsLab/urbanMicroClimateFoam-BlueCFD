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

\*---------------------------------------------------------------------------*/

#include "blendingLayerVelFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::blendingLayerVelFvPatchVectorField::
blendingLayerVelFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    inputTimeStep(3600),
    Utarget_Field(0)
{
}


Foam::blendingLayerVelFvPatchVectorField::
blendingLayerVelFvPatchVectorField
(
    const blendingLayerVelFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    inputTimeStep(ptf.inputTimeStep),
    Utarget_Field(ptf.Utarget_Field)
{}


Foam::blendingLayerVelFvPatchVectorField::
blendingLayerVelFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    inputTimeStep(readLabel(dict.lookup("inputTimeStep"))),
    Utarget_Field(Zero)
{
    fvPatchVectorField::operator=(vectorField("value", dict, p.size()));
}

Foam::blendingLayerVelFvPatchVectorField::
blendingLayerVelFvPatchVectorField
(
    const blendingLayerVelFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    inputTimeStep(ptf.inputTimeStep),
    Utarget_Field(ptf.Utarget_Field)
{}

Foam::blendingLayerVelFvPatchVectorField::
blendingLayerVelFvPatchVectorField
(
    const blendingLayerVelFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    inputTimeStep(ptf.inputTimeStep),
    Utarget_Field(ptf.Utarget_Field)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::blendingLayerVelFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    scalar timeValue = this->db().time().value();
    scalar timeIndex = this->db().time().timeIndex();
    word boundaryName = this->patch().name();
      
    if (timeIndex == 1)
    {
        label moduloTest = int(timeValue/inputTimeStep); 
        word UtargetFile = "Utarget_" + boundaryName + "/Utarget_" + boundaryName +"_" + name(moduloTest*inputTimeStep);
        IOList<vector> Utarget(
            IOobject
            (
                UtargetFile,
                db().time().caseConstant(),
                db(),
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )            
        );
        /////interpolate between two input files if necessary/////
        if (timeValue/inputTimeStep - moduloTest > 0)
        {
            word UtargetFile_B = "Utarget_" + boundaryName + "/Utarget_" + boundaryName +"_" + name(moduloTest*inputTimeStep+inputTimeStep);
            IOList<vector> Utarget_B(
                IOobject
                (
                    UtargetFile_B,
                    db().time().caseConstant(),
                    db(),
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )            
            );
            scalar ratio = (timeValue-moduloTest*inputTimeStep)/inputTimeStep;
            Utarget = Utarget*(1-ratio) + Utarget_B*(ratio);       
        }
        //////////////////////////////////////////////////////////
        
        if (Pstream::parRun()) //in parallel runs, need to read correct values from global wdrFile
        {
            const labelIOList& localFaceProcAddr //global addresses for local faces, includes also internal faces
            (
                IOobject
                (
            	"faceProcAddressing",
            	this->patch().boundaryMesh().mesh().facesInstance(),
            	this->patch().boundaryMesh().mesh().meshSubDir,
            	this->patch().boundaryMesh().mesh(),
            	IOobject::MUST_READ,
            	IOobject::NO_WRITE
                )
            );

            label startFace = this->patch().start(); //local startFace for this patch
            label nFaces = this->patch().size(); //local total face number for this patch

            List<label> globalFaceAddr_;
            globalFaceAddr_.setSize(nFaces);
            forAll(globalFaceAddr_, i) //global address for local faces, only for this patch
            {
                globalFaceAddr_[i] = localFaceProcAddr[startFace + i] - 1; //subtracted 1 to get global index, as localFaceProcAddr starts from 1, not 0
            }

            label minGlobalFaceAddr_ = gMin(globalFaceAddr_); //get the minimum global address for this patch = startFace in global patch

            List<vector> Utarget_;
            Utarget_.setSize(nFaces);
            forAll(Utarget_, i) //read correct values from global wdrFile
            {
                Utarget_[i] = Utarget[localFaceProcAddr[startFace + i] - 1 - minGlobalFaceAddr_]; //subtracted 1 to get global index, as localFaceProcAddr starts from 1, not 0
            }
            
            Utarget_Field = Utarget_;
        }
        else
        {
            Utarget_Field = Utarget;
        }      
    }

    /////ensure mass balance over all lateral boundaries//////    
    word patches[] = {"west", "east", "north", "south"};
    List<scalar> massFlux_WENS(4);
    List<scalar> corrFactor_WENS(4);
    forAll(massFlux_WENS,i)
    {
        label patchId = this->patch().boundaryMesh().findPatchID(patches[i]);
        massFlux_WENS[i] = gSum(-1* this->patch().boundaryMesh().mesh().boundary()[patchId].lookupPatchField<surfaceScalarField, scalar>("phi"));        
    }
    Pstream::listCombineGather(massFlux_WENS, sumOp<scalar>());
    Pstream::listCombineScatter(massFlux_WENS);

    scalar corrFactor;
    forAll(massFlux_WENS,i)
    {
        if (boundaryName == patches[i])
        {
            corrFactor = 1 - sign(massFlux_WENS[i]) * (sum(massFlux_WENS))/(sum(mag(massFlux_WENS)));
        }
    }
    //////////////////////////////////////////////////////////
    
    operator==(Utarget_Field*corrFactor);

    fixedValueFvPatchVectorField::updateCoeffs();
}

void Foam::blendingLayerVelFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("inputTimeStep")
        << inputTimeStep << token::END_STATEMENT << nl;    
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        blendingLayerVelFvPatchVectorField
    );
}

// ************************************************************************* //
