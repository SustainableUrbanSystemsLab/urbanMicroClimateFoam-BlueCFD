Info<< "Reading simulation-physics settings" << endl;
IOdictionary physicsSettings
(
	IOobject
	(
		"physicsSettings",
		runTime.constant(),
		runTime.db(),
		IOobject::MUST_READ,
		IOobject::NO_WRITE
	)
);

const bool solveHeatEquation = physicsSettings.lookupOrDefault("solveHeatEquation", false);	
Info << "solveHeatEquation: " << solveHeatEquation << "\n" << endl;