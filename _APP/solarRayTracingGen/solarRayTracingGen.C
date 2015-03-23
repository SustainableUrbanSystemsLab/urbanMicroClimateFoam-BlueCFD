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

Application
    solarRayTracingGen

Description
    Aytac Kubilay, 2015, Empa
	Based on viewFactorsGen

\*---------------------------------------------------------------------------*/


#include "argList.H"
#include "fvMesh.H"
#include "Time.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "distributedTriSurfaceMesh.H"
#include "cyclicAMIPolyPatch.H"
#include "triSurfaceTools.H"
#include "mapDistribute.H"

#include "OFstream.H"
#include "meshTools.H"
#include "plane.H"
#include "uindirectPrimitivePatch.H"
#include "DynamicField.H"
#include "IFstream.H"
#include "unitConversion.H"

#include "mathematicalConstants.H"
#include "scalarMatrices.H"
#include "CompactListList.H"
#include "labelIOList.H"
#include "labelListIOList.H"
#include "scalarListIOList.H"
#include "scalarIOList.H"
#include "vectorIOList.H"

#include "singleCellFvMesh.H"
#include "IOdictionary.H"
#include "fixedValueFvPatchFields.H"
#include "wallFvPatch.H"

#include "unitConversion.H"

using namespace Foam;

void writeRays
(
    const fileName& fName,
    const pointField& compactCf,
    const pointField& myFc,
    const labelListList& visibleFaceFaces
)
{
    OFstream str(fName);
    label vertI = 0;

    Pout<< "Dumping rays to " << str.name() << endl;

    forAll(myFc, faceI)
    {
        const labelList visFaces = visibleFaceFaces[faceI];
        forAll(visFaces, faceRemote)
        {
            label compactI = visFaces[faceRemote];
            const point& remoteFc = compactCf[compactI];

            meshTools::writeOBJ(str, myFc[faceI]);
            vertI++;
            meshTools::writeOBJ(str, remoteFc);
            vertI++;
            str << "l " << vertI-1 << ' ' << vertI << nl;
        }
    }
    string cmd("objToVTK " + fName + " " + fName.lessExt() + ".vtk");
    Pout<< "cmd:" << cmd << endl;
    system(cmd);
}


scalar calculateViewFactorFij
(
    const vector& i,
    const vector& j,
    const vector& dAi,
    const vector& dAj
)
{
    vector r = i - j;
    scalar rMag = mag(r);
    scalar dAiMag = mag(dAi);
    scalar dAjMag = mag(dAj);

    vector ni = dAi/dAiMag;
    vector nj = dAj/dAjMag;
    scalar cosThetaJ = mag(nj & r)/rMag;
    scalar cosThetaI = mag(ni & r)/rMag;

    return
    (
        (cosThetaI*cosThetaJ*dAjMag*dAiMag)
       /(sqr(rMag)*constant::mathematical::pi)
    );
}


void insertMatrixElements
(
    const globalIndex& globalNumbering,
    const label fromProcI,
    const labelListList& globalFaceFaces,
    const scalarListList& viewFactors,
    scalarSquareMatrix& matrix
)
{
    forAll(viewFactors, faceI)
    {
        const scalarList& vf = viewFactors[faceI];
        const labelList& globalFaces = globalFaceFaces[faceI];

        label globalI = globalNumbering.toGlobal(fromProcI, faceI);
        forAll(globalFaces, i)
        {
            matrix[globalI][globalFaces[i]] = vf[i];
        }
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "addRegionOption.H"
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createNamedMesh.H"

    // Read view factor dictionary
    IOdictionary viewFactorDict
    (
       IOobject
       (
            "viewFactorsDict",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
       )
    );

    const bool writeViewFactors =
        viewFactorDict.lookupOrDefault<bool>("writeViewFactorMatrix", false);

    const bool dumpRays =
        viewFactorDict.lookupOrDefault<bool>("dumpRays", false);
		
	//vector sunPos = viewFactorDict.lookup("sunPosVector");	
	vector skyPos = viewFactorDict.lookup("skyPosVector");

    // Read sunPosVector2 list
    vectorIOList sunPosVector2
    (
       IOobject
       (
            "sunPosVector2",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
       )
    );
    scalarIOList Eb
    (
       IOobject
       (
            "Eb",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
       )
    );    

    const label debug = viewFactorDict.lookupOrDefault<label>("debug", 0);

    volScalarField Qr
    (
        IOobject
        (
            "Qr",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        ),
        mesh
    );

    // Read agglomeration map
    labelListIOList finalAgglom
    (
        IOobject
        (
            "finalAgglom",
            mesh.facesInstance(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    // Create the coarse mesh  using agglomeration
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    if (debug)
    {
        Info << "\nCreating single cell mesh..." << endl;
    }

    singleCellFvMesh coarseMesh
    (
        IOobject
        (
            mesh.name(),
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        finalAgglom
    );


    // Calculate total number of fine and coarse faces
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    label nCoarseFaces = 0;      //total number of coarse faces
	label nCoarseFacesAll = 0;   //Also includes non-wall faces with greyDiffusive boundary
    label nFineFaces = 0;        //total number of fine faces

    const polyBoundaryMesh& patches = mesh.boundaryMesh();
    const polyBoundaryMesh& coarsePatches = coarseMesh.boundaryMesh();

    labelList viewFactorsPatches(patches.size());
	labelList howManyCoarseFacesPerPatch(patches.size());
    const volScalarField::GeometricBoundaryField& Qrb = Qr.boundaryField();

    label count = 0;
	label countAll = 0;
    forAll(Qrb, patchI)
    {
        const polyPatch& pp = patches[patchI];
        const fvPatchScalarField& QrpI = Qrb[patchI];

        //if ((isA<fixedValueFvPatchScalarField>(QrpI)) && (pp.size() > 0))
		if ((isA<wallFvPatch>(mesh.boundary()[patchI])) && (pp.size() > 0))
        {
            viewFactorsPatches[count] = QrpI.patch().index();
            nCoarseFaces += coarsePatches[patchI].size();
			nCoarseFacesAll += coarsePatches[patchI].size();
            nFineFaces += patches[patchI].size();
			count ++;
			
			howManyCoarseFacesPerPatch[countAll] = coarsePatches[patchI].size();          
        }
		else if ((isA<fixedValueFvPatchScalarField>(QrpI)) && (pp.size() > 0))
		{
			nCoarseFacesAll += coarsePatches[patchI].size();
			
			howManyCoarseFacesPerPatch[countAll] = coarsePatches[patchI].size();          
		}
		else 
		{
			howManyCoarseFacesPerPatch[countAll] = 0;          
		}
		countAll ++;
    }
    viewFactorsPatches.resize(count--);
	
	Info << "howManyCoarseFacesPerPatch: " << howManyCoarseFacesPerPatch << endl;
	
    // total number of coarse faces
    label totalNCoarseFaces = nCoarseFaces;

    reduce(totalNCoarseFaces, sumOp<label>());

    if (Pstream::master())
    {
        Info << "\nTotal number of coarse faces: "<< totalNCoarseFaces << endl;
    }

    if (Pstream::master() && debug)
    {
        Pout << "\nView factor patches included in the calculation : "
             << viewFactorsPatches << endl;
    }

    // Collect local Cf and Sf on coarse mesh
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    DynamicList<point> localCoarseCf(nCoarseFaces);
    DynamicList<point> localCoarseSf(nCoarseFaces);

    forAll (viewFactorsPatches, i)
    {
        const label patchID = viewFactorsPatches[i];

        const polyPatch& pp = patches[patchID];
        const labelList& agglom = finalAgglom[patchID];
		
        label nAgglom = max(agglom)+1;
        labelListList coarseToFine(invertOneToMany(nAgglom, agglom));
        const labelList& coarsePatchFace = coarseMesh.patchFaceMap()[patchID];

        const pointField& coarseCf = coarseMesh.Cf().boundaryField()[patchID];
        const pointField& coarseSf = coarseMesh.Sf().boundaryField()[patchID];

        forAll(coarseCf, faceI)
        {
            point cf = coarseCf[faceI];
            const label coarseFaceI = coarsePatchFace[faceI];
            const labelList& fineFaces = coarseToFine[coarseFaceI];
            // Construct single face
            uindirectPrimitivePatch upp
            (
                UIndirectList<face>(pp, fineFaces),
                pp.points()
            );

            List<point> availablePoints
            (
                upp.faceCentres().size()
              + upp.localPoints().size()
            );

            SubList<point>
            (
                availablePoints,
                upp.faceCentres().size()
            ).assign(upp.faceCentres());

            SubList<point>
            (
                availablePoints,
                upp.localPoints().size(),
                upp.faceCentres().size()
            ).assign(upp.localPoints());

            point cfo = cf;
            scalar dist = GREAT;
            forAll(availablePoints, iPoint)
            {
                point cfFine = availablePoints[iPoint];
                if(mag(cfFine-cfo) < dist)
                {
                    dist = mag(cfFine-cfo);
                    cf = cfFine;
                }
            }

            point sf = coarseSf[faceI];
            localCoarseCf.append(cf);
            localCoarseSf.append(sf);
        }
    }

	
    // Set up searching engine for obstacles
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #include "searchingEngine.H"
	

    // Determine rays between coarse face centres
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    DynamicList<label> rayStartFace(nCoarseFaces + 0.01*nCoarseFaces);

    DynamicList<label> rayEndFace(rayStartFace.size());

    globalIndex globalNumbering(nCoarseFaces);
	
	
	// Find the bounding box of the domain
    // ######################################
	point min_(point::zero);
	point max_(point::zero);
	for (label i = 1; i < mesh.points().size(); i++)
	{
		min_ = ::Foam::min(min_, mesh.points()[i]);
		max_ = ::Foam::max(max_, mesh.points()[i]);
	}

	// Find the Solar Ray Start Points within domain
    // ######################################	
	List<point> solarStart(localCoarseCf);
	DynamicField<point> solarEnd(solarStart.size());

    // Number of visible faces from local index
    labelList nVisibleFaceFaces(nCoarseFaces, 0);
    labelListList nVisibleFaceFaces3(sunPosVector2.size());

    forAll(sunPosVector2, vectorId)
    {
        labelList nVisibleFaceFaces2(nCoarseFaces, 0);
        
        vector sunPos = sunPosVector2[vectorId];

    	//List<pointIndexHit> hitInfo(1);
    	forAll(solarStart, pointI)
    	{
    		scalar i1 = 0; scalar i2 = 0; scalar i3 = 0;

    		if (sunPos.x() > 0.0)
            {
                i1 = (max_.x() - solarStart[pointI].x())/sunPos.x();
            } 
            else if (sunPos.x() < 0.0)
            {
                i1 = (min_.x() - solarStart[pointI].x())/sunPos.x();
            } 
            else {i1 = VGREAT;}

    		if (sunPos.y() > 0.0)
            {
                i2 = (max_.y() - solarStart[pointI].y())/sunPos.y();
            } 
            else if (sunPos.y() < 0.0)
            {
                i2 = (min_.y() - solarStart[pointI].y())/sunPos.y();
            }
            else{i2 = VGREAT;}

    		if (sunPos.z() > 0.0)
            {
                i3 = (max_.z() - solarStart[pointI].z())/sunPos.z();
            } 
            else if (sunPos.z() < 0.0)
            {
                i3 = (min_.z() - solarStart[pointI].z())/sunPos.z();
            }
            else{i3 = VGREAT;}

    		scalar i = min(i1, min(i2, i3));
    		point solarEndPoint = i*point(sunPos.x(),sunPos.y(),sunPos.z())+solarStart[pointI];
    		solarEnd.append(solarEndPoint);
    	}
    	
    	// Collect Cf and Sf on coarse mesh
        // #############################################

        List<pointField> remoteCoarseCf_(Pstream::nProcs());
        //List<pointField> remoteCoarseSf(Pstream::nProcs());

        remoteCoarseCf_[Pstream::myProcNo()] = solarEnd;
        //remoteCoarseSf[Pstream::myProcNo()] = localCoarseSf; 
    	
        List<pointField> localCoarseCf_(Pstream::nProcs());
        localCoarseCf_[Pstream::myProcNo()] = solarStart; 

        List<pointField> localCoarseSf_(Pstream::nProcs());
        localCoarseSf_[Pstream::myProcNo()] = localCoarseSf;  	

    	// Collect remote Cf and Sf on fine mesh
        // #############################################

     /*   List<pointField> remoteFineCf(Pstream::nProcs());
        List<pointField> remoteFineSf(Pstream::nProcs());

        remoteCoarseCf[Pstream::myProcNo()] = solarEnd;
        remoteCoarseSf[Pstream::myProcNo()] = localCoarseSf;*/
    	
        // Distribute local coarse Cf and Sf for shooting rays
        // #############################################

        Pstream::gatherList(remoteCoarseCf_);
        Pstream::scatterList(remoteCoarseCf_);
    	Pstream::gatherList(localCoarseCf_);
    	Pstream::scatterList(localCoarseCf_);
        Pstream::gatherList(localCoarseSf_);
        Pstream::scatterList(localCoarseSf_);
        //Pstream::gatherList(remoteCoarseSf);
        //Pstream::scatterList(remoteCoarseSf);


        // Return rayStartFace in local index and rayEndFace in global index
        // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

        #include "shootRays.H"


        forAll(rayStartFace, i)
        {
            nVisibleFaceFaces2[rayStartFace[i]]++;
        }

        labelListList visibleFaceFaces(nCoarseFaces);

        label nViewFactors = 0;
        forAll(nVisibleFaceFaces2, faceI)
        {
            visibleFaceFaces[faceI].setSize(nVisibleFaceFaces2[faceI]);
            nViewFactors += nVisibleFaceFaces2[faceI];
        }

    	//Info << "nVisibleFaceFaces2: " << nVisibleFaceFaces2 << endl;
        nVisibleFaceFaces3[vectorId] = nVisibleFaceFaces2;

        rayStartFace.clear();
        rayEndFace.clear();

    }

    //Info << "nVisibleFaceFaces3: " << nVisibleFaceFaces3 << endl;

    // Fill local view factor matrix
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    labelListIOList sunVisibleOrNot
    (
        IOobject
        (
            "sunVisibleOrNot",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        sunPosVector2.size()
    ); 
    scalarListIOList sunViewCoeff
    (
        IOobject
        (
            "sunViewCoeff",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        sunPosVector2.size()
    );
    scalarIOList skyViewCoeff
    (
        IOobject
        (
            "skyViewCoeff",
            mesh.facesInstance(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        nCoarseFacesAll
    );	
    
    labelList dummy(nCoarseFacesAll, -1);
    scalarList dummy2(nCoarseFacesAll, 0.0);
    forAll(sunVisibleOrNot, vectorId)
    {
        sunVisibleOrNot[vectorId] = dummy;
        sunViewCoeff[vectorId] = dummy2;
	}
	forAll(skyViewCoeff, faceI)
	{
		skyViewCoeff[faceI] = 0;
	}

	scalar cosPhi = 0;
	scalar radAngleBetween = 0;
	scalar degAngleBetween = 0;
	
	label faceNo = 0;
	label i = 0;
	label j = 0;
	label k = 0;

    forAll(sunPosVector2, vectorId)
    {    
        vector sunPos = sunPosVector2[vectorId];

    	forAll(viewFactorsPatches, patchID)
    	{
    		while (i < viewFactorsPatches[patchID])
    		{
    			while (j < howManyCoarseFacesPerPatch[i])
    			{
    				sunVisibleOrNot[vectorId][k] = 0;
    				k++;
    				j++;
    			}
    			j = 0;
    			i++;
    		} 
    		
    		while (j < howManyCoarseFacesPerPatch[i])
    		{
    			sunVisibleOrNot[vectorId][k] = nVisibleFaceFaces3[vectorId][faceNo];
    			
    			cosPhi = (localCoarseSf[faceNo] & sunPos)/(mag(localCoarseSf[faceNo])*mag(sunPos) + SMALL);
    			sunViewCoeff[vectorId][k] = nVisibleFaceFaces3[vectorId][faceNo]*mag(cosPhi)*Eb[vectorId];

    			cosPhi = (localCoarseSf[faceNo] & skyPos)/(mag(localCoarseSf[faceNo])*mag(skyPos) + SMALL);
    			radAngleBetween = Foam::acos( min(max(cosPhi, -1), 1) );
    			degAngleBetween = radToDeg(radAngleBetween);
    			if (degAngleBetween == 180){degAngleBetween=0;}
    			else if (degAngleBetween > 90){degAngleBetween-=90;}
    			skyViewCoeff[k] = 1-0.5*(degAngleBetween/90);			
    			
    			k++;
    			j++;
    			faceNo++;
    		}
    	}
        i = 0;
        j = 0;
        k = 0;
        faceNo = 0;
    }
	Info << "sunVisibleOrNot: " << sunVisibleOrNot << endl;	

	Info << "localCoarseCf: " << localCoarseCf << endl;	
	Info << "localCoarseSf: " << localCoarseSf << endl;
	
	Info << "sunViewCoeff: " << sunViewCoeff << endl;	
	Info << "skyViewCoeff: " << skyViewCoeff << endl;	

	sunVisibleOrNot.write();
	sunViewCoeff.write();	
	skyViewCoeff.write();	

    Info<< "End\n" << endl;
    return 0;
}


// ************************************************************************* //
