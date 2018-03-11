//update_sequence = 1; 
/*pc.dimensions().reset(dimless);
Krel.dimensions().reset(dimless);
K_v.dimensions().reset(dimless);
K_pt.dimensions().reset(dimless);
Crel.dimensions().reset(dimless);
//ws.dimensions().reset(dimless);
//#include "updatebuildingMaterials.H" 
//ws.dimensions().reset(dimensionSet( 1, -3, 0, 0, 0, 0, 0));
Krel.dimensions().reset(dimensionSet( 0, 0, 1, 0, 0, 0, 0));    
K_v.dimensions().reset(dimensionSet( 0, 0, 1, 0, 0, 0, 0));  
K_pt.dimensions().reset(dimensionSet( 1, -1, -1, -1, 0, 0, 0));            
pc.dimensions().reset(dimensionSet( 1, -1, -2, 0, 0, 0, 0));    
Crel.dimensions().reset(dimensionSet( 0, -2, 2, 0, 0, 0, 0)); */

//pc.correctBoundaryConditions();
//ws.correctBoundaryConditions();

//Crel.correctBoundaryConditions();
/*fvScalarMatrix dummy(fvm::laplacian(Crel,pc,"laplacian(Krel,pc)"));
Info << "000: " << dummy.boundaryCoeffs() << endl;
fvScalarMatrix test(fvm::laplacian(Crel,pc,"laplacian(Krel,pc)"));
Info << "AAA: " << test.boundaryCoeffs() << endl;*/

Krel.boundaryField().evaluate();
K_v.correctBoundaryConditions();
Crel.correctBoundaryConditions();

/*fvScalarMatrix test2(fvm::laplacian(Crel,pc,"laplacian(Krel,pc)"));
Info << "BBB: " << test2.boundaryCoeffs() << endl;*/

//K_pt.correctBoundaryConditions();//pc.write(); Krel.write(); K_v.write(); w.write(); 

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
//pcEqn.relax();
pcEqn.solve(); 

