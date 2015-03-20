simpleControlFluid simpleControlFluid(fluidRegions[i]);

label fluidIteration = 0;

while (simpleControlFluid.loop())
{
//  Pressure-velocity SIMPLE corrector

    p_rgh.storePrevIter();
    rho.storePrevIter();
    {
        #include "UEqn.H"
        #include "EEqn.H"
        #include "pEqn.H"
    }

    turb.correct();

    fluidIteration++;
    Info << nl;
}

Info<< nl << fluidRegions[i].name() << " region solution converged in "
    << fluidIteration << " iterations" << nl << endl;