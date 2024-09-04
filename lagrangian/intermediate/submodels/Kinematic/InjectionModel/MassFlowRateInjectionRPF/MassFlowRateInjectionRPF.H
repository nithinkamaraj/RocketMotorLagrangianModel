/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
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

Class
    Foam::MassFlowRateInjectionRPF

Group
    grpLagrangianIntermediateInjectionSubModels

Description
    Patch injection, by using patch flow rate to determine massLoadingRatio and
    velocity.

    User specifies:
      - Total mass to inject
      - Name of patch
      - Injection duration
      - Injection target massLoadingRatio/carrier volume flow rate

    Properties:
      - Initial parcel velocity given by local flow velocity
      - Parcel diameters obtained by distribution model
      - Parcels injected randomly across the patch

SourceFiles
    MassFlowRateInjectionRPF.C

\*---------------------------------------------------------------------------*/

#ifndef MassFlowRateInjectionRPF_H
#define MassFlowRateInjectionRPF_H

#include "InjectionModel.H"
#include "patchInjectionBaseRPF.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

class distributionModel;

/*---------------------------------------------------------------------------*\
                   Class MassFlowRateInjectionRPF Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class MassFlowRateInjectionRPF
:
    public InjectionModel<CloudType>,
    public patchInjectionBaseRPF
{
    // Private data

        //- Name of carrier (mass or volume) flux field
        const word phiName_;

        //- Name of carrier density field
        const word rhoName_;

        //- Injection duration [s]
        scalar duration_;

        //- massLoadingRatio profile of particle volume to carrier volume [-]
        autoPtr<Function1<scalar>> massLoadingRatio_;

        //- Parcels to introduce per unit volume flow rate m3 [n/m3]
        const scalar nParticle_;

        //- Parcel size distribution model
        const autoPtr<distributionModel> sizeDistribution_;


public:

    //- Runtime type information
    TypeName("MassFlowRateInjectionRPF");


    // Constructors

        //- Construct from dictionary
        MassFlowRateInjectionRPF
        (
            const dictionary& dict,
            CloudType& owner,
            const word& modelName
        );

        //- Construct copy
        MassFlowRateInjectionRPF(const MassFlowRateInjectionRPF<CloudType>& im);

        //- Construct and return a clone
        virtual autoPtr<InjectionModel<CloudType>> clone() const
        {
            return autoPtr<InjectionModel<CloudType>>
            (
                new MassFlowRateInjectionRPF<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~MassFlowRateInjectionRPF();


    // Member Functions

        //- Inherit updateMesh from patchInjectionBaseRPF
        using patchInjectionBaseRPF::updateMesh;

        //- Set injector locations when mesh is updated
        virtual void updateMesh();

        //- Return the end-of-injection time
        scalar timeEnd() const;

        // //- Return the total volumetric flow rate across the patch [m3/s]
        // virtual scalar flowRate() const;

        //- Return the total mass flow rate across the patch [kg/s]
        virtual scalar massFlowRate() const;

        //- Number of parcels to introduce relative to SOI
        virtual label parcelsToInject(const scalar time0, const scalar time1);

        //- Volume of parcels to introduce relative to SOI
        virtual scalar volumeToInject(const scalar time0, const scalar time1);


        // Injection geometry

            //- Inherit setPositionAndCell from patchInjectionBaseRPF
            using patchInjectionBaseRPF::setPositionAndCell;

            //- Set the injection position and owner cell, tetFace and tetPt
            virtual void setPositionAndCell
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                vector& position,
                label& cellOwner,
                label& tetFacei,
                label& tetPti
            );

            virtual void setProperties
            (
                const label parcelI,
                const label nParcels,
                const scalar time,
                typename CloudType::parcelType& parcel
            );

            //- Flag to identify whether model fully describes the parcel
            virtual bool fullyDescribed() const;

            //- Return flag to identify whether or not injection of parcelI is
            //  permitted
            virtual bool validInjection(const label parcelI);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "MassFlowRateInjectionRPF.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
