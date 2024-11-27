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
@Inbook{Kashani1986,
author="Kashani, A.
and Van Sciver, S. W.",
editor="Fast, R. W.",
title="Steady State Forced Convection Heat Transfer in He II",
bookTitle="Advances in Cryogenic Engineering: Volume 31",
year="1986",
publisher="Springer US",
address="Boston, MA",
pages="489--498",
abstract="A study of forced convection heat transfer in superfluid helium (He II) is initiated to better understand the physical behavior of this process and to compare it with the more familiar He II heat transfer mechanism of internal convection. An experimental assembly is designed to achieve fluid flow by a motor-driven hydraulic pump which utilizes two stainless steel bellows. Each bellows is connected to one end of a copper tube, 3 mm in diameter and 2 m long. The system allows measurements of one dimensional heat and mass transfer where the measured quantities include: temperature profile and pressure drop. The variable quantities are the helium bath temperature, flow velocity and heat input. The helium bath is held at 1.8 K and under saturation pressure. The flow tube is heated at the middle and the flow velocity is varied up to 97 cm/s. The helium pressure is monitored at both ends of the tube and a friction factor is estimated for He II. Temperature measurements are made at seven evenly spaced locations along the tube. The experimental temperature profile is compared with a numerical solution of an analytical model developed for the problem under study.",
isbn="978-1-4613-2213-9",
doi="10.1007/978-1-4613-2213-9_56",
url="https://doi.org/10.1007/978-1-4613-2213-9_56"
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
    #include "volHeatingSource.H"

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
            );

            TEqn.solve();
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
