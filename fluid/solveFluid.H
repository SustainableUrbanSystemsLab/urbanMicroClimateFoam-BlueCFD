simpleControlFluid simpleControlFluid(fluidRegions[i]);

TimeState pts(runTime); //store time state

label fluidIteration = 1;

while ((simpleControlFluid.run(runTime) || fluidIteration <= minFluidIteration) && fluidIteration <= maxFluidIteration)
{
//  Pressure-velocity SIMPLE corrector
    Info << "Internal iteration number: " << fluidIteration << endl;

    gra.correct(thermo.T(), w, U, Sh, Sw);

    if(!simple.lookupOrDefault<bool>("frozenFlow", 0))
    {
        p_rgh.storePrevIter();
        rho.storePrevIter();
        {
            #include "UEqn.H"
            #include "EEqn.H"
            #include "pEqn.H"
            #include "wEqn.H"
        }
    turb.correct();
    }

    Info << nl; 
    fluidIteration++;
}

Info<< nl << fluidRegions[i].name() << " region solution converged in "
    << fluidIteration-1 << " iterations" << nl << endl;

Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
rad.correct();
Info << "Updating short-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
sol.correct();

runTime.TimeState::operator=(pts); //restore time state

