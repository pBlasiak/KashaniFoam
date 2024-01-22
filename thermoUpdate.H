laminarTransport.correct();
rho = laminarTransport.rhoHe();
cp = laminarTransport.cpHe();

gradT = fvc::grad(T);	
magGradT = max(mag(gradT), dimensionedScalar("small", dimTemperature/dimLength, SMALL)); 
magGradT2 = pow(magGradT,2);

laminarTransport.calckHe();
kHe = laminarTransport.kHe();
