        dimensionedScalar muair("muair",dimensionSet(1,-1,-1,0,0,0,0),scalar(1.8e-5));
        dimensionedScalar Pr("Pr",dimensionSet(0,0,0,0,0,0,0),scalar(0.7));
        dimensionedScalar cp("cp",dimensionSet(0,2,-2,1,0,0,0),scalar(1005));

        surfaceScalarField heatFlux
        (
		 fvc::interpolate(turb.alphaEff())*cp*fvc::snGrad(thermo.T())
        );

        const surfaceScalarField::GeometricBoundaryField& patchHeatFlux=
           heatFlux.boundaryField();

        Info<< "\nmappedWall heat fluxes" << endl;
        forAll(patchHeatFlux, patchi)
        {
            {
                Info<< mesh.boundary()[patchi].name()
                    << " "
                    << gSum
                       (
                           mesh.magSf().boundaryField()[patchi]
                          *patchHeatFlux[patchi]
                       )
                    << endl;
            }
        }
        Info<< endl;

        forAll(wallHeatFlux.boundaryField(), patchi)
        {
            wallHeatFlux.boundaryField()[patchi] = patchHeatFlux[patchi];
        }
