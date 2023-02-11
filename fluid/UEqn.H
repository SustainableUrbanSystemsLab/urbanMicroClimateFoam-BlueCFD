    // Solve the Momentum equation

    MRF.correctBoundaryVelocity(U);

    tmp<fvVectorMatrix> tUEqn
    (
        fvm::div(phi, U)
      + MRF.DDt(rho, U)
      + turb.divDevRhoReff(U)
     ==
      - fvm::SuSp(gra.Cf()*rho*mag(U), U)
      + fvOptions(rho, U)
      - fvm::SuSp(veg.Cf()*rho*mag(U), U) // vegetation momentum source
    );
    fvVectorMatrix& UEqn = tUEqn.ref();

    if (blending)
    {
        UEqn -= rho*bL.bL_USource(U);
    }

    UEqn.relax();

    fvOptions.constrain(UEqn);

    solve
    (
        UEqn
     ==
        fvc::reconstruct
        (
            (
              - ghf*fvc::snGrad(rho)
              - fvc::snGrad(p_rgh)
            )*mesh.magSf()
        )
    );

    fvOptions.correct(U);
