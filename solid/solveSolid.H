{
    TimeState pts(runTime); //store time state

    runTime.setTime(pts.timeOutputValue()-pts.deltaTValue(),pts.timeIndex()-1);
    runTime.setDeltaT(pts.deltaTValue() * initialSolidTimestepFactor);

    scalar solidInternalTime = 0;
    scalar timeToOutput = pts.deltaTValue();
    bool timeStepDecrease = false;

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
            if (gMin(Ts) <= 0 || gMin(Ts.boundaryField()) <= 0)
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
            if (timeToOutput >= 0.0)
            {
                runTime.setTime(runTime.value()+runTime.deltaT().value(),runTime.timeIndex()+1);
                //runTime++;  //using this creates problems when writeInterval != 1
            }
            #include "setSolidRegionDeltaT.H"
        }

        Info << "timeToOutput: " << timeToOutput << endl;
    }

    runTime.TimeState::operator=(pts); //restore time state
}


