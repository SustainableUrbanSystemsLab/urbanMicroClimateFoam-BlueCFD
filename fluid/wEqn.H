{
    //dimensionedScalar Dm("Dm",dimensionSet(0,2,-1,0,0,0,0),scalar(0.282e-4)); 
    //scalar Sct = 0.7;

    fvScalarMatrix wEqn
    (
        fvm::ddt(rho,w)
      + fvm::div(phi, w)
//      - fvm::laplacian(turb.alphaEff(), w) //In reactingParcelFimFoam, they use alphaEff, but in others Mass Diffusivity = muEff/Sct      
      - fvm::laplacian(turb.mut()/0.7, w)
//      - fvm::laplacian(rho*Dm, w)      
    );
    wEqn.relax();

    wEqn.solve();

}
