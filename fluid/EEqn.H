{

    volScalarField& he = thermo.he();

    fvScalarMatrix EEqn
    (
        fvm::div(phi, he)
      + (
            he.name() == "e"
          ? fvc::div(phi, volScalarField("Ekp", 0.5*magSqr(U) + p/rho))
          : fvc::div(phi, volScalarField("K", 0.5*magSqr(U)))
        )
      + fluidThermophys.divq(he)
     ==
        rho*(U&g)
      + rad.Sh(thermo, he)
      + sol.Sh(thermo, he)
      + gra.Sh()
      + veg.Sh() // vegetation energy source term
      + fvOptions(rho, he)
    );

    if (blending)
    {
        EEqn -= rho*thermo.Cp()*bL.bL_TSource(thermo.T());
    }

    EEqn.relax();

    fvOptions.constrain(EEqn);

    EEqn.solve();

    fvOptions.correct(he);

    thermo.correct();

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
