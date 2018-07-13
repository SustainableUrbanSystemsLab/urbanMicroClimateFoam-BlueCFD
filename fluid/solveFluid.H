simpleControlFluid simpleControlFluid(fluidRegions[i]);

const scalar storeFluidDeltaT = runTime.deltaT().value();
const scalar storeFluidTimeStep = runTime.value();
const label storeFluidTimeIndex = runTime.timeIndex();
const dimensionedScalar storeFluidEndTime = runTime.endTime();
runTime.setTime(0,0);
runTime.setDeltaT(1);

label fluidIteration = 1;
while ((simpleControlFluid.loop(runTime) || fluidIteration <= minFluidIteration) && fluidIteration <= maxFluidIteration)
{
//  Pressure-velocity SIMPLE corrector
    Info << "Internal iteration number: " << fluidIteration << endl;

    p_rgh.storePrevIter();
    rho.storePrevIter();
    {
        #include "UEqn.H"
        #include "EEqn.H"
        #include "pEqn.H"
        #include "wEqn.H"
    }

    turb.correct();

    Info << nl; 
    fluidIteration++;
}

Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
rad.correct();
Info << "Updating short-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
sol.correct();

runTime.setDeltaT(storeFluidDeltaT);
runTime.setTime(storeFluidTimeStep,storeFluidTimeIndex);
runTime.setEndTime(storeFluidEndTime);
