/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionary.H"
#include "Ostream.H"
#include "Pstream.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 31 "/home/openfoam/run/bL/0/U/boundaryField/inlet/#codeStream"
#include "fvCFD.H"
//}}} end codeInclude

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    void codeStream_bc1feaa5ab31bcf953a1e5c8d0f433f397efaa55
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 48 "/home/openfoam/run/bL/0/U/boundaryField/inlet/#codeStream"
const IOdictionary& d = static_cast<const IOdictionary&>
                (
                    dict.parent().parent()
                );

                const fvMesh& mesh = refCast<const fvMesh>(d.db());
                const label id = mesh.boundary().findPatchID("inlet");
                const fvPatch& patch = mesh.boundary()[id];

                // Initialize the vector field
                vectorField U(patch.size(), vector(0, 0, 0));

                scalar U_0;
                const scalar U_inf = 14.03838;

                const scalar a1 =  4.409e+23;
                const scalar a2 = -4.73e+21;
                const scalar a3 =  2.053e+19;
                const scalar a4 = -4.542e+16;
                const scalar a5 =  5.125e+13;
                const scalar a6 = -1.981e+10;
                const scalar a7 = -1.552e+07;
                const scalar a8 =  2.048e+04;
                const scalar a9 =  -0.02949;

                const scalar b1 = -5.964e+12;
                const scalar b2 =  1.425e+12;
                const scalar b3 = -1.447e+11;
                const scalar b4 =  8.13e+09;
                const scalar b5 = -2.764e+08;
                const scalar b6 =  5.839e+06;
                const scalar b7 = -7.656e+04;
                const scalar b8 =  711.4;
                const scalar b9 =  7.028;

                forAll(U, i)
                {
                    const scalar y = patch.Cf()[i][1];

                    if (y <= 0.002534939443643) {
                        U_0 = a1*pow(y,8) + a2*pow(y,7) + a3*pow(y,6) + a4*pow(y,5) + \
                              a5*pow(y,4) + a6*pow(y,3) + a7*pow(y,2) + a8*y + a9;
                        U[i] = vector(U_0, 0., 0.);
                    } else if ( y > 0.002534939443643 && y <= 0.050054661656196) {
                        U_0 = b1*pow(y,8) + b2*pow(y,7) + b3*pow(y,6) + b4*pow(y,5) + \
                              b5*pow(y,4) + b6*pow(y,3) + b7*pow(y,2) + b8*y + b9;
                        U[i] = vector(U_0, 0., 0.);
                    } else {
                        U[i] = vector(U_inf, 0., 0.);
                    }
                }

                writeEntry(os, "", U);
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

