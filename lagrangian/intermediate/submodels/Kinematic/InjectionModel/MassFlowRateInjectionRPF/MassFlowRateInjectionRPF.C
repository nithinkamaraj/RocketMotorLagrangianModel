/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2015-2020 OpenCFD Ltd.
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

#include "MassFlowRateInjectionRPF.H"
#include "distributionModel.H"
#include "mathematicalConstants.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MassFlowRateInjectionRPF<CloudType>::MassFlowRateInjectionRPF
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName,typeName),
    patchInjectionBaseRPF(owner.mesh(), this->coeffDict().getWord("patch")),
    phiName_(this->coeffDict().template getOrDefault<word>("phi", "phi")),
    rhoName_(this->coeffDict().template getOrDefault<word>("rho", "rho")),
    duration_(this->coeffDict().getScalar("duration")),
    massLoadingRatio_
    (
        Function1<scalar>::New
        (
            "massLoadingRatio",
            this->coeffDict(),
            &owner.mesh()
        )
    ),
    nParticle_
    (
        this->coeffDict().getScalar("nParticle")
    ),
    sizeDistribution_
    (
        distributionModel::New
        (
            this->coeffDict().subDict("sizeDistribution"),
            owner.rndGen()
        )
    )
{
    // Convert from user time to reduce the number of time conversion calls
    const Time& time = owner.db().time();
    duration_ = time.userTimeToTime(duration_);
    massLoadingRatio_->userTimeToTime(time);

    patchInjectionBaseRPF::updateMesh(owner.mesh());

    // Re-initialise total mass/volume to inject to zero
    // - will be reset during each injection
    this->volumeTotal_ = 0.0;
    this->massTotal_ = 0.0;
}


template<class CloudType>
Foam::MassFlowRateInjectionRPF<CloudType>::MassFlowRateInjectionRPF
(
    const MassFlowRateInjectionRPF<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    patchInjectionBaseRPF(im),
    phiName_(im.phiName_),
    rhoName_(im.rhoName_),
    duration_(im.duration_),
    massLoadingRatio_(im.massLoadingRatio_.clone()),
    nParticle_(im.nParticle_),
    sizeDistribution_(im.sizeDistribution_.clone())
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::MassFlowRateInjectionRPF<CloudType>::~MassFlowRateInjectionRPF()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::MassFlowRateInjectionRPF<CloudType>::updateMesh()
{
    patchInjectionBaseRPF::updateMesh(this->owner().mesh());
}


template<class CloudType>
Foam::scalar Foam::MassFlowRateInjectionRPF<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}

/// gives the massFlow rate
template<class CloudType>
Foam::scalar Foam::MassFlowRateInjectionRPF<CloudType>::massFlowRate() const
{
   const polyMesh& mesh = this->owner().mesh();

   const surfaceScalarField& phi =
        mesh.lookupObject<surfaceScalarField>(phiName_);

   const scalarField& phip = phi.boundaryField()[patchId_];

   scalar massFlowRateIn = 0.0;

   const volScalarField& rho =
          mesh.lookupObject<volScalarField>(rhoName_);

   const scalarField& rhop = rho.boundaryField()[patchId_];

   massFlowRateIn = max(0.0, -sum(phip));


  reduce(massFlowRateIn, sumOp<scalar>());

       Info<< "massFlowRate   "<< massFlowRateIn << endl;
    return massFlowRateIn;
}


template<class CloudType>
Foam::label Foam::MassFlowRateInjectionRPF<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))

    {
        scalar dt = time1 - time0;


        scalar c = massLoadingRatio_->value(0.5*(time0 + time1));



        scalar d = sizeDistribution_->sample();


        scalar massp = nParticle_*this->owner().constProps().rho0()*pi/6.0*pow3(d);



        scalar nParcels = c*massFlowRate()*dt/massp;





        Info<< "time0                               "<< time0 << endl;
        Info<< "time1                               "<< time1 << endl;
        Info<< "dt                                  "<< dt << endl;
        Info<< "massLoadingRatio                        "<< c << endl;
        Info<< "diameter                            "<< d << endl;
        Info<< "no of a parcel                      "<< nParticle_ << endl;
        Info<< "no of parcels in this time step     "<< nParcels << endl;



        Random& rnd = this->owner().rndGen();

        label nParcelsToInject = floor(nParcels);

        // Inject an additional parcel with a probability based on the
        // remainder after the floor function
        if
        (
            nParcelsToInject > 0
         && (
               nParcels - scalar(nParcelsToInject)
             > rnd.globalPosition(scalar(0), scalar(1))
            )
        )
        {
            ++nParcelsToInject;
        }


        Info<<" no of parcels injected      "<< nParcelsToInject << endl;




        return nParcelsToInject;
    }

    return 0;
}


template<class CloudType>
Foam::scalar Foam::MassFlowRateInjectionRPF<CloudType>::volumeToInject
(
    const scalar time0,
    const scalar time1
)
{
    scalar mass = 0.0;
    scalar volume = 0.0;

    if ((time0 >= 0.0) && (time0 < duration_))
    {
        scalar c = massLoadingRatio_->value(0.5*(time0 + time1));

        mass = c*(time1 - time0)*massFlowRate();
    }

    this->massTotal_ = mass;

    this->volumeTotal_ = mass/this->owner().constProps().rho0();

     volume = this->volumeTotal_;
    return volume;
}


template<class CloudType>
void Foam::MassFlowRateInjectionRPF<CloudType>::setPositionAndCell
(
    const label,
    const label,
    const scalar,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{

    patchInjectionBaseRPF::setPositionAndCell
    (
        this->owner().mesh(),
        this->owner().rndGen(),
        position,
        cellOwner,
        tetFacei,
        tetPti
    );
}


template<class CloudType>
void Foam::MassFlowRateInjectionRPF<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // Set particle velocity to carrier velocity
    parcel.U() = vector(0,0,0); //this->owner().U()[parcel.cell();

    // Set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::MassFlowRateInjectionRPF<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::MassFlowRateInjectionRPF<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
