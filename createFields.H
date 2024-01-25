Info<< "Reading field T\n" << endl;

volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField T
(
    IOobject
    (
        "T",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"

singlePhaseHeliumTransportModel laminarTransport(U, phi);

laminarTransport.correct();

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	laminarTransport.rhoHe()
);

volScalarField rhon
(
    IOobject
    (
        "rhon",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	rho
);
	rhon = laminarTransport.rhon();

volScalarField rhos
(
    IOobject
    (
        "rhos",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	rho
);
	rhos = rho - rhon;

volScalarField sHe
(
    IOobject
    (
        "s",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	laminarTransport.sHe()
);

volScalarField AGM
(
    IOobject
    (
        "AGM",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	laminarTransport.AGMHe()
);

volScalarField cp
(
    IOobject
    (
        "cp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	laminarTransport.cpHe()
);

volScalarField onebyf
(
    IOobject
    (
        "onebyf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	laminarTransport.onebyf()
);

volScalarField kHe
(
    IOobject
    (
        "kHe",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
    dimensionedScalar("small", dimensionSet(1,1,-3,-1,0,0,0), SMALL)
);

volVectorField gradT
(
    IOobject
    (
        "gradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("small", dimensionSet(0,-1,0,1,0,0,0), vector(0,0,0))
);
gradT = fvc::grad(T);

volScalarField magGradT
(
    IOobject
    (
        "magGradT",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
    dimensionedScalar("small", dimensionSet(0,-1,0,1,0,0,0), SMALL)
);
magGradT = mag(gradT);

volScalarField magGradT2
(
    IOobject
    (
        "magGradT2",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
	mesh,
    dimensionedScalar("small", dimensionSet(0,-2,0,2,0,0,0), SMALL)
);


#include "createFvOptions.H"
