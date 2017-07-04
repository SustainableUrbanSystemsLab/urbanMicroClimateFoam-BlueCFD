simpleControlFluid simpleControlFluid(fluidRegions[i]);

label fluidIteration = 0;

while (simpleControlFluid.loop() && fluidIteration <= maxFluidIteration)
{
//  Pressure-velocity SIMPLE corrector
    fluidIteration++;
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
}

Info<< nl << fluidRegions[i].name() << " region solution converged in "
    << fluidIteration << " iterations" << nl << endl;

Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
rad.correct();
Info << "Updating short-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
sol.correct();

//#include "wallMoistureFlux.H"
//#include "wallHeatFlux.H"
