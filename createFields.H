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

//- Helium thermal conductivity models enumerations
//  kHe = (1/(f|gradT|^2))^(1/3)
enum class HeliumThermalConductivityModel
{
	Hepak,
	Kitamura,
	Suekane
};

static const Enum<HeliumThermalConductivityModel> HeliumThermalConductivityModelNames_
{
    { HeliumThermalConductivityModel::Hepak, "Hepak" },
    { HeliumThermalConductivityModel::Kitamura, "Kitamura" },
    { HeliumThermalConductivityModel::Suekane, "Suekane" },
};

HeliumThermalConductivityModel kHeModel
(
	HeliumThermalConductivityModelNames_.get("kHeModel", laminarTransport)
);

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

volScalarField Qsource
(
    IOobject
    (
        "Qsource",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
	mesh,
    dimensionedScalar("small", dimensionSet(1,-1,-3,0,0,0,0), 0)
);

dimensionedScalar q("q", dimPower/dimLength/dimLength, laminarTransport);
dimensionedScalar d("pipeDiameter", dimLength, laminarTransport);
dimensionedScalar heaterLocation("heaterLocation", dimLength, laminarTransport);
dimensionedScalar Q = constant::mathematical::pi*d*d/4*q;
Info<< "\nApplied heat transfer rate Q = " << Q << endl;

forAll(Qsource, celli)
{
	scalar x = Qsource.mesh().C()[celli].component(vector::X);
	if (x > heaterLocation.value())
	{
		Qsource[celli] = Q.value()/Qsource.mesh().V()[celli];
		Info<< "Qsource = " << Qsource[celli] << ", V = " << Qsource.mesh().V()[celli] << endl;
		break;
	}
}

//Info<< "Source term vector: " << Qsource << endl;

//volScalarField DT
//(
//    IOobject
//    (
//        "DT",
//        runTime.timeName(),
//        mesh,
//        IOobject::READ_IF_PRESENT,
//        IOobject::AUTO_WRITE
//    ),
//    mesh,
//    dimensionedScalar(dimViscosity, Zero)
//);
//
//if (!DT.headerOk())
//{
//    IOdictionary transportProperties
//    (
//        IOobject
//        (
//            "transportProperties",
//            runTime.constant(),
//            mesh,
//            IOobject::MUST_READ_IF_MODIFIED,
//            IOobject::NO_WRITE
//        )
//    );
//    DT = dimensionedScalar("DT", dimViscosity, transportProperties);
//}

#include "createFvOptions.H"
