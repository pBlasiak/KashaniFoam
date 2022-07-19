/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    KashaniFoam

Group
    grpBasicSolvers

Description
    Solver for one-dimensional forced convection of He II.


    \heading Solver details
    The solver solves 1D energy equation in He II.  The
    equation is given by:
TODO: poprawic rownanie

    \f[
        \ddt{T}  = \div \left( D_T \grad T \right)
    \f]

    Where:
    \vartable
        T     | Scalar field which is solved for, e.g. temperature
        D_T   | Diffusion coefficient
    \endvartable

    \heading Required fields
    \plaintable
        T     | Scalar field which is solved for, e.g. temperature
    \endplaintable

References
@article{doi:10.1080/10407788908944714,
author = { A.   Kashani  and  S. W.   Van Sciver  and  J. C.   Strikwerda },
title = {NUMERICAL SOLUTION OF FORCED CONVECTION HEAT TRANSFER IN He II},
journal = {Numerical Heat Transfer, Part A: Applications},
volume = {16},
number = {2},
pages = {213-228},
year  = {1989},
publisher = {Taylor & Francis},
doi = {10.1080/10407788908944714},

URL = { 
        https://doi.org/10.1080/10407788908944714
    
},
eprint = { 
        https://doi.org/10.1080/10407788908944714
}
}

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "singlePhaseHeliumTransportModel.H"
#include "fvOptions.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "1D energy equation for forced convection in He II."
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    simpleControl simple(mesh);

    #include "createFields.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating temperature distribution\n" << endl;

    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
            #include "thermoUpdate.H"

            fvScalarMatrix TEqn
            (
                (rho*U*cp & fvc::grad(T))   // rhoPhiCp??
              - fvm::laplacian(kHe, T) 
             ==
			    Qsource
		//	  + fvOptions(T)
            );

         //   fvOptions.constrain(TEqn);
            TEqn.solve();
          //  fvOptions.correct(T);
        }

		Info<< "min/max(T) = " << min(T).value() << ", "
    	    << max(T).value() <<endl;

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
