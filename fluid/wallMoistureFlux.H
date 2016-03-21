    dimensionedScalar Dm("Dm",dimensionSet(0,2,-1,0,0,0,0),scalar(2.5e-5)); 
    scalar Sct = 0.7;
    dimensionedScalar rhoair("Dm",dimensionSet(1,-3,0,0,0,0,0),scalar(1.2)); 

        surfaceScalarField moistureFlux
        (
            (rhoair*Dm + fvc::interpolate(turb.mut())/Sct)*fvc::snGrad(w)
        );

        const surfaceScalarField::GeometricBoundaryField& patchMoistureFlux =
           moistureFlux.boundaryField();

        Info<< "\nmappedWall moisture fluxes" << endl;
        forAll(patchMoistureFlux, patchi)
        {
            {
                Info<< mesh.boundary()[patchi].name()
                    << " "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchMoistureFlux[patchi]
                       )
                    << endl;
            }
        }
        Info<< endl;

        forAll(wallMoistureFlux.boundaryField(), patchi)
        {
            wallMoistureFlux.boundaryField()[patchi] = patchMoistureFlux[patchi];
        }
