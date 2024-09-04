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

#include "CellMassFlowRateInjection.H"
#include "distributionModel.H"
#include "mathematicalConstants.H"
#include "surfaceFields.H"


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CellMassFlowRateInjection<CloudType>::CellMassFlowRateInjection
(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName
)
:
    InjectionModel<CloudType>(dict, owner, modelName,typeName),
    CellInjectionBase(owner.mesh(), this->coeffDict().getWord("patch")),
    phiName_(this->coeffDict().template getOrDefault<word>("phi.gas", "phi.gas")),
    rhoName_(this->coeffDict().template getOrDefault<word>("rho.gas", "rho.gas")),
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
    ),
    time0_(),
    time1_()
{
     // Convert from user time to reduce the number of time conversion calls
     const Time& time = owner.db().time();
     duration_ = time.userTimeToTime(duration_);

     this->volumeTotal_ = 0.0;
     this->massTotal_ = 0.0;

     updateMesh();

/////////////////////////////////////////////////////////
     const polyMesh& mesh = this->owner().mesh();
     const polyPatch& patch = mesh.boundaryMesh()[patchId_];
     const pointField& points = patch.points();
     DynamicList<face> tris(5);
     forAll(patch, facei) // what to do when you restart ? (Revisit when problem arrises!)
     {
          const face& f = patch[facei];
          tris.clear();
          f.triangles(points, tris);

          forAll(tris, i)
          {
            trinPr_.append(0.0);
          }
     }


     Info<<" processor no        "<<Pstream::myProcNo() <<endl;
     ///Info<<"size of remPr       "<<  remPr[Pstream::myProcNo()].size()<<endl;

   // listlistPr_[Pstream::myProcNo()].transfer(trinPi_);

   Info<<"   Reminder parcel list is intialized  "<< endl;
/////////////////////////////////////////////////////////////
}


template<class CloudType>
Foam::CellMassFlowRateInjection<CloudType>::CellMassFlowRateInjection
(
    const CellMassFlowRateInjection<CloudType>& im
)
:
    InjectionModel<CloudType>(im),
    CellInjectionBase(im),
    phiName_(im.phiName_),
    rhoName_(im.rhoName_),
    duration_(im.duration_),
    massLoadingRatio_(im.massLoadingRatio_.clone()),
    nParticle_(im.nParticle_),
    sizeDistribution_(im.sizeDistribution_.clone()),
    time0_(im.time0_),
    time1_(im.time1_)
{}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CellMassFlowRateInjection<CloudType>::~CellMassFlowRateInjection()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CellMassFlowRateInjection<CloudType>::updateMesh()
{
  Info<<"////////////////////////////// UPDATE MESH ////////////////////////////////"<<endl;

  const polyMesh& mesh = this->owner().mesh();
  const polyPatch& patch = mesh.boundaryMesh()[patchId_];
  const pointField& points = patch.points();

  const volVectorField& U =
       mesh.lookupObject<volVectorField>("U.gas");
  const vectorField& Up = U.boundaryField()[patchId_];

//Info <<"U.gas"<<min(Up)<<" - "<<max(Up)<<endl;

  const volScalarField& rho =
       mesh.lookupObject<volScalarField>("rho.gas");
  const scalarField& rhop = rho.boundaryField()[patchId_];

  // Info <<"rho.gas"<<min(rhop)<<" - "<<max(rhop)<<endl;

  const scalar time = this->owner().db().time().value();
  Info << " Time      " << time << endl;

  const scalar t0 = this->timeStep0_ - this->SOI_;
  const scalar t1 = time - this->SOI_;

  scalar nParticle    =     nParticle_;                                         // number of particle in a parcel
  scalar d            =     sizeDistribution_->sample();                        // dimaeter of particle
  scalar mlr          =     massLoadingRatio_->value(0.5*(t0 + t1));            // Returns mass loading ratio - constant
  scalar rho0         =     this->owner().constProps().rho0();                  // particle density
  scalar massp        =     rho0*pi*pow3(d)/6;                                  // mass of one particle
  scalar dt           =     t1 - t0 ;                                           // time step

  cellOwners_ = patch.faceCells();
  Info<< " size of the patch cells   " << cellOwners_.size() << endl;

  DynamicList<label> triToFace(2*patch.size()); // face id associated with the triangle
  DynamicList<scalar> trinPi(2*patch.size());
  DynamicList<scalar> triMF(2*patch.size());
  DynamicList<face> triFace(2*patch.size());  // traingle id in a assciated face
  DynamicList<face> tris(5);

  forAll(patch, facei)
  {
    const face& f = patch[facei];

    tris.clear();
    f.triangles(points, tris);

    const vector& Uface = Up[facei];

    forAll(tris, i)
    {
        triToFace.append(facei);
        triFace.append(tris[i]);

        // number of parcels entering the triangle in the present time step
        scalar nPi =  (mlr * tris[i].mag(points) * mag(Uface.x()) * rhop[facei] *dt)/(massp*nParticle);
        trinPi.append(nPi);

        // Mass flow rate per traingle
        triMF.append(tris[i].mag(points) * mag(Uface.x()) * rhop[facei]);
    }

  }


Info<<" No of processers used   "<< Pstream::nProcs() <<endl;
Info<< "size of list            "<<triFace.size()<<endl;

  // calculating number of parcels to be injected and reminder parcels
  DynamicList<scalar> trinP(2*patch.size());
  DynamicList<scalar> trinPr(2*patch.size());

  for(label i = 0; i < triFace.size(); i++)
  {
      // total number of particles to be injected (includes previous time step particles)
      scalar nP = trinPi[i];

      if (trinPr_.size() > 0)
      {
          nP += trinPr_[i];
      }
      // scalar nP = trinPi[i] + listlistPr_[Pstream::myProcNo()][i];

      // Flooring the number of particles to inject
      trinP.append(floor(nP));

      // Finding reminder parcel count and storing
      scalar nPr = nP - floor(nP);
      trinPr.append(nPr);
  }

  triFace_.transfer(triFace);
  triToFace_.transfer(triToFace);
  trinP_.transfer(trinP);
  triMF_.transfer(triMF);
  trinPr_.transfer(trinPr);

  Info<<"trinP                    "<< sum(trinP_)<<endl;
  Info<<"trinPr                   "<< sum(trinPr_)<<endl;

Info<<"/////////////////  UPDATE MESH COMPLETED //////////////////"<<endl;

  // listlistPr_[Pstream::myProcNo()].transfer(trinPr);
////////////////////////////////////////////////////

    // calculating patch normal and patch area
    const scalarField magSf(mag(patch.faceAreas()));
    patchArea_ = sum(magSf);
    patchNormal_ = patch.faceAreas()/magSf;
    reduce(patchArea_, sumOp<scalar>());
}


template<class CloudType>
Foam::scalar Foam::CellMassFlowRateInjection<CloudType>::timeEnd() const
{
    return this->SOI_ + duration_;
}

/// gives the massFlow rate
template<class CloudType>
Foam::scalar Foam::CellMassFlowRateInjection<CloudType>::massFlowRate() const
{
   const polyMesh& mesh           =   this->owner().mesh();
   const surfaceScalarField& phi  =   mesh.lookupObject<surfaceScalarField>(phiName_);
   const scalarField& phip        =   phi.boundaryField()[patchId_];
   scalar massFlowRateIn = 0.0;

     const volScalarField& rho =
          mesh.lookupObject<volScalarField>("rho.gas");
     const scalarField& rhop = rho.boundaryField()[patchId_];

   massFlowRateIn = max(0.0, -sum(rhop*phip));
   reduce(massFlowRateIn, sumOp<scalar>());

   scalar triFaceMassFlowRate =0.0;
   triFaceMassFlowRate = max(0.0, sum(this->triMF_));
   reduce(triFaceMassFlowRate, sumOp<scalar>());


  Info << "massFlowRate             " << massFlowRateIn << endl;
  Info << "triFaceMassFlowRate         "<< triFaceMassFlowRate << endl;

   return  triFaceMassFlowRate  ;
}


template<class CloudType>
Foam::label Foam::CellMassFlowRateInjection<CloudType>::parcelsToInject
(
    const scalar time0,
    const scalar time1
)
{
    if ((time0 >= 0.0) && (time0 < duration_))
    {
        scalar nParcelsToInject = 0.0;
        nParcelsToInject = max(0.0, sum(this->trinP_));
        reduce(nParcelsToInject, sumOp<scalar>());


Info << " no of parcels injected :    "<< nParcelsToInject << endl;



        return nParcelsToInject ;
    }

    return 0 ;
}

template<class CloudType>
Foam::scalar Foam::CellMassFlowRateInjection<CloudType>::volumeToInject
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
void Foam::CellMassFlowRateInjection<CloudType>::setPositionAndCell
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

    CellInjectionBase::setPositionAndCell
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
void Foam::CellMassFlowRateInjection<CloudType>::setProperties
(
    const label,
    const label,
    const scalar,
    typename CloudType::parcelType& parcel
)
{
    // Set particle velocity to carrier velocity
    parcel.U() = vector(0,0,0);             //this->owner().U()[parcel.cell()];

    // Set particle diameter
    parcel.d() = sizeDistribution_->sample();
}


template<class CloudType>
bool Foam::CellMassFlowRateInjection<CloudType>::fullyDescribed() const
{
    return false;
}


template<class CloudType>
bool Foam::CellMassFlowRateInjection<CloudType>::validInjection(const label)
{
    return true;
}


// ************************************************************************* //
