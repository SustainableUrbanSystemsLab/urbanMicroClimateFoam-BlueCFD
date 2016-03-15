
   Info<< "Reading field w\n" << endl;
    volScalarField w
    (
        IOobject
        (
            "w",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );
	
    Info<< "Reading field mut\n" << endl;

    volScalarField mut
    (
        IOobject
        (
            "mut",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );		
