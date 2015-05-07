{
    instantList timeDirs = timeSelector::select0(runTime, args);
    runTime.setTime(timeDirs.last(),timeDirs.size()-1);

    runTime.setDeltaT(storeFluidDeltaT * initialSolidTimestepFactor);

    scalar solidInternalTime = runTime.deltaT().value();

    while ( solidInternalTime <= storeFluidDeltaT )
    {
        volScalarField T_old = thermo.T();
        bool timeStepDecrease = false;

        runTime++; 
        Info << nl << "Time = " << runTime.timeName() << ", deltaT = " << runTime.deltaT().value() << endl;  
        
        for (int nInternalIter=1; nInternalIter<=nInternalIterMax; nInternalIter++)
        {
            for (int nonOrth=0; nonOrth<=nNonOrthCorr; nonOrth++)
            {
                fvScalarMatrix hEqn
                (
                    fvm::ddt(betav*rho, h)
                  - fvm::laplacian(betav*alpha, h, "laplacian(alpha,h)")
                 ==
                    fvOptions(rho, h)
                );

                hEqn.relax();

                fvOptions.constrain(hEqn);

                hEqn.solve();

                fvOptions.correct(h);
            }

            thermo.correct();
            Info<< "Min/max T:" << min(thermo.T()) << ' ' << max(thermo.T()) << endl;

            //convergence test
            scalar maxChangeT = gMax(mag(thermo.T().internalField()-T_old.internalField())); 
            if(maxChangeT < 1e-1 && nInternalIter>=2) //force at least 1 internal iteration
            {
                T_old = thermo.T();
                Info << "Total internal iterations: " << nInternalIter << endl;
                timeStepDecrease = false;
                break;
            }
            else if (nInternalIter == nInternalIterMax) //nonlinear iteration reached maximum
            {
                T_old = thermo.T();
                Info
                   << "Internal iteration didn't converge !"
                   << endl
                   << "maxChangeT: " << maxChangeT << endl
                   << "Total internal iterations: " << nInternalIter << endl;
                timeStepDecrease = true;
            }   
            else //not converged nor reached the maximum iteration yet, continue
            {
                T_old = thermo.T();
                Info << "maxChangeT: " << maxChangeT << endl;
            }               

            //update thermal radiation fluxes
            forAll(fluidRegions, i)
            {
                Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name();
                radiation::radiationModel& rad = radiation[i];
                rad.correct();
            }
        }

        scalar timeToOutput = storeFluidDeltaT - solidInternalTime;
        if (timeToOutput > 0.0)
        {
            #include "readSolidTimeControls.H"
            #include "setSolidRegionDeltaT.H"
        }

        solidInternalTime += runTime.deltaT().value();
    }

    runTime.setDeltaT(storeFluidDeltaT);
}


