simpleControlFluid simpleControlFluid(fluidRegions[i]);

label fluidIteration = 1;

while (simpleControlFluid.loop())
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

    fluidIteration++;
    Info << nl;
}

Info<< nl << fluidRegions[i].name() << " region solution converged in "
    << fluidIteration << " iterations" << nl << endl;

Info << "Updating long-wave radiation heat transfer for region: " << fluidRegions[i].name();
rad.correct();
Info << "Updating short-wave radiation heat transfer for region: " << fluidRegions[i].name();
sol.correct();