/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2021 OpenCFD Ltd.
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

#include "SphereDragForce.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class CloudType>
Foam::scalar Foam::SphereDragForce<CloudType>::CdRe(const scalar Re) const
{
    // // (AOB:Eq. 35)
    // if (Re > 1000.0)
    // {
    //     return 0.424*Re;
    // }
    //
    // return 24.0*(1.0 + (1.0/6.0)*pow(Re, 2.0/3.0));

label a = 24.0*(1.0 + (0.15)*pow(Re, 0.687)+(0.0175)*pow(Re, 2.16)/(pow(Re, 1.16)+42500));
 

//  label a = 0.46*Re + 28*pow(Re,0.15);

    return a;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::SphereDragForce<CloudType>::SphereDragForce
(
    CloudType& owner,
    const fvMesh& mesh,
    const dictionary& dict
)
:
    ParticleForce<CloudType>(owner, mesh, dict, typeName, false)
{}


template<class CloudType>
Foam::SphereDragForce<CloudType>::SphereDragForce
(
    const SphereDragForce<CloudType>& df
)
:
    ParticleForce<CloudType>(df)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
Foam::forceSuSp Foam::SphereDragForce<CloudType>::calcCoupled
(
    const typename CloudType::parcelType& p,
    const typename CloudType::parcelType::trackingData& td,
    const scalar dt,
    const scalar mass,
    const scalar Re,
    const scalar muc
) const
{
//     // (AOB:Eq. 34)
//       // label re = td.rhoc()*mag(p.U()-td.Uc())*p.d()/muc;
//         Info<< "Check: mass of parcel -"<< mass<<endl;
//         Info<< "Check: CdRe -     "<<CdRe(Re)<<endl;
//         Info<< "Check: Re -       "<<Re<<endl;
//       //   Info<< "Check: re -"<<re<<endl;
//       //   Info<< "Check: muc -"<<muc<<endl;
//       //   Info<< "Check: parcel density  -"<<p.rho()<<endl;
//         Info<< "Sp calcualtion -  "<<mass*0.75*muc*CdRe(Re)/(p.rho()*sqr(p.d()))<<endl;
// // 

  return forceSuSp(Zero, mass*0.75*muc*CdRe(Re)/(p.rho()*sqr(p.d())));
}


// ************************************************************************* //
