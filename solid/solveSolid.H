{
    TimeState pts(runTime); //store time state

    runTime.setTime(pts.timeOutputValue()-pts.deltaTValue(),pts.timeIndex()-1);
    runTime.setDeltaT(pts.deltaTValue() * initialSolidTimestepFactor);

    scalar solidInternalTime = 0;
    scalar timeToOutput = pts.deltaTValue();
    bool timeStepDecrease = false;
    
    scalar timeAfterLastRadUpdate = 0;
    
    //store old T//////////
    PtrList<volScalarField::Boundary> T_old(fluidRegions.size());
    forAll(fluidRegions, i)
    {
        rhoThermo& thermo = thermoFluid[i];
        T_old.set
        (
            i,
            new volScalarField::Boundary
            (
                thermo.T().boundaryField()
            )
        );
    }
    ///////////////////////

    while ( solidInternalTime < pts.deltaTValue() )
    {         
        Info << nl << "Time = " << runTime.timeName() << ", deltaT = " << runTime.deltaT().value() << endl;  
        Info << "solidInternalTime: " << solidInternalTime << endl;
        //pc.storePrevIter();
        //Ts.storePrevIter();

        #include "updatebuildingMaterials.H"

        //store values from previous timestep (for mixed form moisture equation)
        volScalarField ws_old = ws; //ws_old.boundaryField().updateCoeffs();
        volScalarField pc_old = pc; //pc_old.boundaryField().updateCoeffs();
        volScalarField Ts_old = Ts; //Ts_old.boundaryField().updateCoeffs();

        //store values from previous Picard iteration
        volScalarField ws_n = ws; //ws_n.boundaryField().updateCoeffs();
        volScalarField pc_n = pc; //pc_n.boundaryField().updateCoeffs();
        volScalarField Ts_n = Ts; //Ts_n.boundaryField().updateCoeffs();

        for (int nInternalIter=1; nInternalIter<=nInternalIterMax; nInternalIter++) //starting Picard iteration
        {

            //Moisture transfer////////////
            #include "pcEqn.H"    

            //Firstly, test if all pc values are valid
            if (gMax(pc) >= 0 || gMax(pc.boundaryField()) >= 0)
            {
                Info << "This is going to crash (pc)! Decreasing timestep and reverting fields..." << endl;
                Info << "Error: gMax(pc): " << gMax(pc) << ", gMax(pc.boundaryField()): " << gMax(pc.boundaryField()) << endl;
                timeStepDecrease = true;
                #include "setSolidRegionDeltaT.H"
                #include "revertValues.H"    
                break;
            }
            pc.correctBoundaryConditions();
            ///////////////////////////////

            //Heat transfer////////////////
            #include "TsEqn.H" 
            //Firstly, test if all Ts values are valid
            if (gMin(Ts) <= 0 || gMin(Ts.boundaryField()) <= 0
                || gMax(Ts) >= 423.15 || gMax(Ts.boundaryField()) >= 423.15)
            {           
                Info << "This is going to crash (Ts)! Decreasing timestep and reverting fields..." << endl;
                timeStepDecrease = true;
                #include "setSolidRegionDeltaT.H"
                #include "revertValues.H"  
                break;                      
            }
            Ts.correctBoundaryConditions();
            ///////////////////////////////
            

            //Convergence test/////////////

            //update values for convergence test
            #include "updatebuildingMaterials.H"
          
            //convergence test
            #include "checkConvergence.H"
            ///////////////////////////////
        }                

        if (timeStepDecrease == false) 
        {
            #include "solidContinuityErrs.H"

            solidInternalTime += runTime.deltaT().value();
            timeToOutput = pts.deltaTValue() - solidInternalTime;
            timeAfterLastRadUpdate += runTime.deltaT().value();
            if (timeToOutput >= 0.0)
            {
                runTime.setTime(runTime.value()+runTime.deltaT().value(),runTime.timeIndex()+1);
                //runTime++;  //using this creates problems when writeInterval != 1
            }
            #include "setSolidRegionDeltaT.H"
        }
        
        if (timeToOutput > 0.0 && timeAfterLastRadUpdate >= 600.0)
        {
            //update thermal radiation fluxes
            forAll(fluidRegions, i)
            {
                Info << "Updating T boundary fields..." << endl;
                rhoThermo& thermo = thermoFluid[i];
                thermo.T().correctBoundaryConditions();
            }
            if (vegRegions.size() > 0)
            {
                forAll(vegRegions, i)
                {
                    Info << "Updating T boundary fields..." << endl;
                    volScalarField& vegT = TVeg[i];                   
                    vegT.correctBoundaryConditions();
                    Info << "Updating long-wave radiation heat transfer for region: " << vegRegions[i].name() << endl;
                    radiationModel& rad = radiation2[i];
                    rad.correct();
                }
            }
            else
            {
                forAll(fluidRegions, i)
                {
                    Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
                    radiationModel& rad = radiation[i];
                    rad.correct();
                }
            }
            timeAfterLastRadUpdate = 0;
        }        

        Info << "timeToOutput: " << timeToOutput << endl;
    }

    //restore old T////////
    forAll(fluidRegions, i)
    {
        rhoThermo& thermo = thermoFluid[i];
        
        volScalarField::Boundary& TBf = thermo.T().boundaryFieldRef();
        const volScalarField::Boundary& refTBf = T_old[i];
        
        forAll(TBf, patchi)
        {
            //restore patch value fields
            forAll(TBf[patchi], facei)
            {
                TBf[patchi][facei] = refTBf[patchi][facei];
            }
            
            //restore patch refValue fields
            if
            (   
                isA<mixedFvPatchScalarField>(TBf[patchi])
             && isA<mixedFvPatchScalarField>(refTBf[patchi])
            )
            {
                refCast<mixedFvPatchScalarField>
                (TBf[patchi]).refValue() =
                refCast<const mixedFvPatchScalarField>
                (refTBf[patchi]).refValue();
            }
        }
    }
    forAll(vegRegions, i)
    {
        volScalarField& vegT = TVeg[i];                   
        vegT.correctBoundaryConditions();
    }
    ///////////////////////

    runTime.TimeState::operator=(pts); //restore time state
}


