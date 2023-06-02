laminarTransport.correct();
rho = laminarTransport.rhoHe();
cp = laminarTransport.cpHe();
sHe = laminarTransport.sHe();

gradT = gradP/rho/sHe + fvc::grad(T);	
magGradT = max(mag(gradT), dimensionedScalar("small", dimTemperature/dimLength, SMALL)); 
magGradT2 = pow(magGradT,2);


switch (kHeModel)
{
	case HeliumThermalConductivityModel::Hepak:
		 {
			kHe = pow(laminarTransport.onebyf()/magGradT2, 1./3);  // HEPAK
		 	break;
		 }

	case HeliumThermalConductivityModel::Kitamura:
		 {
			kHe = pow(1./(AGM*rhon/pow(rhos*sHe*T,3)/sHe)/magGradT2, 1./3); // KITAMURA
		 	break;
		 }

	case HeliumThermalConductivityModel::Suekane:
		 {
			const dimensionedScalar sLambda("sLambda", dimensionSet(0,2,-2,-1,0,0,0), 1560);
		 	const dimensionedScalar Ac("Ac", dimensionSet(-1,1,1,0,0,0,0), 1297);
		 	const dimensionedScalar TLambda("TLambda", dimTemperature, 2.168);
		 	const dimensionedScalar c("c", dimTemperature, 0.01);
		 	kHe = pow(1./(Ac/rho/rho/pow(sLambda,4)/pow(TLambda,3)*
		 				pow(pow((T+c)/(TLambda+c),5.7)*(1-pow((T+c)/(TLambda+c),5.7)),-3))
		 			/magGradT2, 1./3); // SUEKANE
		 	//Info<< "f^-1 = " << 1./(Ac/rho/rho/pow(sLambda,4)/pow(TLambda,3)*
		 	//			pow(pow((T+c)/(TLambda+c),5.7)*(1-pow((T+c)/(TLambda+c),5.7)),-3));
		 	break;
		 }

	default:
		 FatalErrorInFunction
    	     << "Unknown helium thermal conductivity model. Valid models are: "
    	     << HeliumThermalConductivityModelNames_ << nl << exit(FatalError);
		 break;
}
