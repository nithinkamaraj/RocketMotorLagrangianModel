Info<< "\nConstructing reacting cloud" << endl;
basicReactingTypeCloud parcels
(
    "reactingCloud1",
    rho,
    U,
    g,
    slgThermo
);



// Mass fraction Feild of parcels
volScalarField thetaP
(
    IOobject
    (
        "thetaP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);


// Volume fraction feild of parcels

volScalarField alphaP
(
    IOobject
    (
        "alphaP",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);


volScalarField SourceCheck
(
    IOobject
    (
        "SourceCheck",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

