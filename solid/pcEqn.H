//Crel.correctBoundaryConditions();
/*fvScalarMatrix dummy(fvm::laplacian(Crel,pc,"laplacian(Krel,pc)"));
Info << "000: " << dummy.boundaryCoeffs() << endl;
fvScalarMatrix test(fvm::laplacian(Crel,pc,"laplacian(Krel,pc)"));
Info << "AAA: " << test.boundaryCoeffs() << endl;*/

fvScalarMatrix pcEqn
(
    Crel*fvm::ddt(pc)                
//    fvm::ddt(Crel,pc)
//    -pc*fvc::grad(Crel)
//    -fvm::Sp(fvm::ddt(Crel),pc)
    ==
    fvm::laplacian(Krel+K_v,pc,"laplacian(Krel,pc)")
    +fvc::laplacian(K_pt,Ts,"laplacian(Krel,pc)")
);
//Info << "AAA: " << pcEqn.boundaryCoeffs() << endl;
pcEqn.solve(); 

