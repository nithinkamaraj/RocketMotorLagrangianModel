/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
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

#include "CellInjectionBase.H"
#include "polyMesh.H"
#include "SubField.H"
#include "Random.H"
#include "triPointRef.H"
#include "volFields.H"
#include "polyMeshTetDecomposition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::CellInjectionBase::CellInjectionBase
(
    const polyMesh& mesh,
    const word& patchName
)
:
    patchName_(patchName),
    patchId_(mesh.boundaryMesh().findPatchID(patchName_)),
    patchArea_(0.0),
    patchNormal_(),
    cellOwners_(),
    triFace_(),
    triToFace_(),
    trinP_(),
    trinPr_(),
    trinPi_(),
    triMF_(),
    listlistPr_()

{
    if (patchId_ < 0)
    {
        FatalErrorInFunction
            << "Requested patch " << patchName_ << " not found" << nl
            << "Available patches are: " << mesh.boundaryMesh().names() << nl
            << exit(FatalError);
    }


}


Foam::CellInjectionBase::CellInjectionBase(const CellInjectionBase& pib)
:
    patchName_(pib.patchName_),
    patchId_(pib.patchId_),
    patchArea_(pib.patchArea_),
    patchNormal_(pib.patchNormal_),
    cellOwners_(pib.cellOwners_),
    triFace_(pib.triFace_),
    triToFace_(pib.triToFace_),
    trinP_(pib.trinP_),
    trinPr_(pib.trinPr_),
    trinPi_(pib.trinPi_),
    triMF_(pib.triMF_),
    listlistPr_(pib.listlistPr_)



{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::CellInjectionBase::updateMesh(const polyMesh& mesh)
{
  // Check Cellmass flow rate injection
}


Foam::label Foam::CellInjectionBase::setPositionAndCell
(
    const fvMesh& mesh,
    const scalar fraction01,
    Random& rnd,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti

)
{
    label facei = -1;
    label trii = 0;



// Pout<<" size of parcle List "<< trinP_.size()<< endl;

    if (cellOwners_.size() > 0)
    {
          // Find corresponding decomposed face triangle
          for(label i = 0; i < trinP_.size(); i++)
          {
              if ( trinP_[i] > 0 )
              {
                 trii = i;
                 trinP_[i] = trinP_[i] - 1;

          // Pout<<" processor check"<< endl;


                 break;
              }
          }



    //      Pout<< "///////////////// Cell postion set //////////////////"<<endl;

          //  Info<< "Parcel location:- triface index " << trii<<endl;
          // Set cellOwner
          facei = triToFace_[trii];
          cellOwner = cellOwners_[facei];

          // Find random point in triangle
          const polyPatch& patch = mesh.boundaryMesh()[patchId_];
          const pointField& points = patch.points();
          const face& tf = triFace_[trii];
          const triPointRef tri(points[tf[0]], points[tf[1]], points[tf[2]]);
          const point pf(tri.randomPoint(rnd));

          // Position perturbed away from face (into domain)
          const scalar a = rnd.position(scalar(0.1), scalar(0.5));
          const vector& pc = mesh.cellCentres()[cellOwner];
          const vector d =
              mag((pf - pc) & patchNormal_[facei])*patchNormal_[facei];

          position = pf - a*d;

          // Try to find tetFacei and tetPti in the current position
          mesh.findTetFacePt(cellOwner, position, tetFacei, tetPti);

          // tetFacei and tetPti not found, check if the cell has changed
          if (tetFacei == -1 ||tetPti == -1)
          {
              mesh.findCellFacePt(position, cellOwner, tetFacei, tetPti);
          }

          // Both searches failed, choose a random position within
          // the original cell
          if (tetFacei == -1 ||tetPti == -1)
          {
              // Reset cellOwner
              cellOwner = cellOwners_[facei];
              const scalarField& V = mesh.V();

              // Construct cell tet indices
              const List<tetIndices> cellTetIs =
                  polyMeshTetDecomposition::cellTetIndices(mesh, cellOwner);

              // Construct cell tet volume fractions
              scalarList cTetVFrac(cellTetIs.size(), Zero);
              for (label teti=1; teti<cellTetIs.size()-1; teti++)
              {
                  cTetVFrac[teti] =
                      cTetVFrac[teti-1]
                    + cellTetIs[teti].tet(mesh).mag()/V[cellOwner];
              }
              cTetVFrac.last() = 1;

              // Set new particle position
              const scalar volFrac = rnd.sample01<scalar>();
              label teti = 0;
              forAll(cTetVFrac, vfI)
              {
                  if (cTetVFrac[vfI] > volFrac)
                  {
                      teti = vfI;
                      break;
                  }
              }
              position = cellTetIs[teti].tet(mesh).randomPoint(rnd);
              tetFacei = cellTetIs[teti].face();
              tetPti = cellTetIs[teti].tetPt();
          }
       }
    else
    {
        cellOwner = -1;
        tetFacei = -1;
        tetPti = -1;

        // Dummy position
        position = pTraits<vector>::max;
    }
  //
  // }
  //   else
  //   {
  //       cellOwner = -1;
  //       tetFacei = -1;
  //       tetPti = -1;
  //
  //       // Dummy position
  //       position = pTraits<vector>::max;
  //   }
  //

    return facei;
}


Foam::label Foam::CellInjectionBase::setPositionAndCell
(
    const fvMesh& mesh,
    Random& rnd,
    vector& position,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
  scalar fraction01 = rnd.globalSample01<scalar>();
    return setPositionAndCell
    (
        mesh,
        fraction01,
        rnd,
        position,
        cellOwner,
        tetFacei,
        tetPti
    );
}





// Foam::label Foam::CellInjectionBase::whichProc(const scalar fraction01) const
// {
//     const scalar areaFraction = fraction01*(sum(triMF_));
//
//     // Determine which processor to inject from
//     forAllReverse(triProcP_, i)
//     {
//         if (areaFraction >= triProcP_[i])
//         {
//             return i;
//
//
//           //  Info <<" processor i    "<< i << endl;
//         }
//     }
//
//     return 0;
// }
//



// ************************************************************************* //
