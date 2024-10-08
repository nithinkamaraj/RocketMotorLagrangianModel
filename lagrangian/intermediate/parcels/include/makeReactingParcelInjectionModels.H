/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#ifndef makeReactingParcelInjectionModels_H
#define makeReactingParcelInjectionModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CellZoneInjection.H"
#include "ConeInjection.H"
#include "ConeNozzleInjection.H"
#include "FieldActivatedInjection.H"
#include "ManualInjection.H"
#include "NoInjection.H"
#include "PatchInjection.H"
#include "PatchFlowRateInjection.H"
#include "MassFlowRateInjection.H"
#include "MassFlowRateInjectionRPF.H"
#include "ReactingLookupTableInjection.H"
#include "CellMassFlowRateInjection.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingParcelInjectionModels(CloudType)                           \
                                                                               \
    makeInjectionModel(CloudType);                                             \
    makeInjectionModelType(CellZoneInjection, CloudType);                      \
    makeInjectionModelType(ConeInjection, CloudType);                          \
    makeInjectionModelType(ConeNozzleInjection, CloudType);                    \
    makeInjectionModelType(FieldActivatedInjection, CloudType);                \
    makeInjectionModelType(ManualInjection, CloudType);                        \
    makeInjectionModelType(NoInjection, CloudType);                            \
    makeInjectionModelType(PatchInjection, CloudType);                         \
    makeInjectionModelType(PatchFlowRateInjection, CloudType);                 \
    makeInjectionModelType(MassFlowRateInjection, CloudType);                  \
    makeInjectionModelType(ReactingLookupTableInjection, CloudType);           \
    makeInjectionModelType(MassFlowRateInjectionRPF, CloudType);               \
    makeInjectionModelType(CellMassFlowRateInjection, CloudType);



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
