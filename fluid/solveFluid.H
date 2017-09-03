simpleControlFluid simpleControlFluid(fluidRegions[i]);

simplifiedVegetationModel vegetation(U, thermo.T(), w, Tl); // Vegetation model

label fluidIteration = 1;

while ((simpleControlFluid.loop() || fluidIteration <= minFluidIteration) && fluidIteration <= maxFluidIteration)
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

    // Solve vegetation energy balance
    vegetation.solve(U, thermo.T(), w, Tl);

    turb.correct();

    Info << nl;
    fluidIteration++;
}

Info<< nl << fluidRegions[i].name() << " region solution converged in "
    << fluidIteration-1 << " iterations" << nl << endl;

Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
rad.correct();
Info << "Updating short-wave radiation heat transfer for region: " << fluidRegions[i].name() << endl;
sol.correct();
