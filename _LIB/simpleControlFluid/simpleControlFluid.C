/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2018 OpenFOAM Foundation
     \\/     M anipulation  |
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

\*---------------------------------------------------------------------------*/

#include "simpleControlFluid.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(simpleControlFluid, 0);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::simpleControlFluid::simpleControlFluid(fvMesh& mesh, const word& algorithmName)
:
    fluidSolutionControl(mesh, algorithmName),
    singleRegionConvergenceControl
    (
        static_cast<singleRegionSolutionControl&>(*this)
    )
{
    read();
    printResidualControls();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::simpleControlFluid::~simpleControlFluid()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

bool Foam::simpleControlFluid::read()
{
    return fluidSolutionControl::read() && readResidualControls();
}


bool Foam::simpleControlFluid::run(Time& time)
{
    read();

    if (converged())
    {
        return false;
    }
    else
    {
        storePrevIterFields();
        time.setDeltaT(0);
        time++; //necessary for iter().stream() to get the correct iteration for convergence control - ayk
        return true;
    }
}

bool Foam::simpleControlFluid::converged()
{
    if
    (
        control_.time().timeIndex() != control_.time().startTimeIndex()
     && criteriaSatisfied()
    )
    {
//        Info<< nl << control_.algorithmName() << " solution converged in "
//            << control_.time().timeName() << " iterations" << nl << endl;
        return true;
    }

    return false;
}

// ************************************************************************* //
