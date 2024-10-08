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

#include "patchInjectionBaseRPF.H"
#include "polyMesh.H"
#include "SubField.H"
#include "Random.H"
#include "triPointRef.H"
#include "volFields.H"
#include "polyMeshTetDecomposition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::patchInjectionBaseRPF::patchInjectionBaseRPF
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
    triCumulativeMagSf_(),
    sumTriMagSf_(Pstream::nProcs() + 1, Zero)
{
    if (patchId_ < 0)
    {
        FatalErrorInFunction
            << "Requested patch " << patchName_ << " not found" << nl
            << "Available patches are: " << mesh.boundaryMesh().names() << nl
            << exit(FatalError);
    }

    updateMesh(mesh);
}


Foam::patchInjectionBaseRPF::patchInjectionBaseRPF(const patchInjectionBaseRPF& pib)
:
    patchName_(pib.patchName_),
    patchId_(pib.patchId_),
    patchArea_(pib.patchArea_),
    patchNormal_(pib.patchNormal_),
    cellOwners_(pib.cellOwners_),
    triFace_(pib.triFace_),
    triToFace_(pib.triToFace_),
    triCumulativeMagSf_(pib.triCumulativeMagSf_),
    sumTriMagSf_(pib.sumTriMagSf_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::patchInjectionBaseRPF::updateMesh(const polyMesh& mesh)
{


Info<<"//////////////// UPDATE MESH ////////////////"<<endl;



    // Set/cache the injector cells
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];
    const pointField& points = patch.points();

        const volVectorField& U =
             mesh.lookupObject<volVectorField>("U");
        const vectorField& Up = U.boundaryField()[patchId_];

        const volScalarField& rho =
             mesh.lookupObject<volScalarField>("rho");
        const scalarField& rhop = rho.boundaryField()[patchId_];



    cellOwners_ = patch.faceCells();

    // Triangulate the patch faces and create addressing
    DynamicList<label> triToFace(2*patch.size());
    DynamicList<scalar> triMagSf(2*patch.size());
    DynamicList<face> triFace(2*patch.size());
    DynamicList<face> tris(5);

    // Set zero value at the start of the tri area list
    triMagSf.append(0.0);

    forAll(patch, facei)
    {
        const face& f = patch[facei];

        tris.clear();
        f.triangles(points, tris);

        forAll(tris, i)
        {
            triToFace.append(facei);
            triFace.append(tris[i]);
            triMagSf.append(tris[i].mag(points)*mag(Up[facei])*rhop[facei]);
        }
    }

Info<< "size of list  "<<triFace.size()<<endl;


    sumTriMagSf_ = Zero;
    sumTriMagSf_[Pstream::myProcNo() + 1] = sum(triMagSf);


    Pstream::listCombineAllGather(sumTriMagSf_, maxEqOp<scalar>());

//    Pout<<"sum of proc area"<<sumTriMagSf_<<endl;




        for (label i = 1; i < triMagSf.size(); i++)
        {
            triMagSf[i] += triMagSf[i-1];
        }

    // Transfer to persistent storage
    triFace_.transfer(triFace);
    triToFace_.transfer(triToFace);
    triCumulativeMagSf_.transfer(triMagSf);


  //  Pout<< "  triCumulativeMagSf_        " << tri<<endl;


    // Convert sumTriMagSf_ into cumulative sum of areas per proc
//Info<< "sumTriMagSf   "<< sumTriMagSf_<< endl;


    for (label i = 1; i < sumTriMagSf_.size(); i++)
    {
        sumTriMagSf_[i] += sumTriMagSf_[i-1];
    }


// Pout<< "  sumTriMagSf_        " << sumTriMagSf_ <<endl;

   const scalarField magUp(mag(Up));
    const scalarField magSf(mag(patch.faceAreas()));
    patchArea_ = sum(magSf*magUp*rhop);
    patchNormal_ = patch.faceAreas()/magSf;
    reduce(patchArea_, sumOp<scalar>());
}


Foam::label Foam::patchInjectionBaseRPF::setPositionAndCell
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

    if (cellOwners_.size() > 0)
    {
        // Determine which processor to inject from
        const label proci = whichProc(fraction01);




        if (Pstream::myProcNo() == proci)
        {
            const scalar areaFraction = fraction01*patchArea_;

            // Find corresponding decomposed face triangle
            label trii = 0;
            scalar offset = sumTriMagSf_[proci];
            forAllReverse(triCumulativeMagSf_, i)
            {
                if (areaFraction > triCumulativeMagSf_[i] + offset)
                {
                    trii = i;


                    break;
                }
            }



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
    }
    else
    {
        cellOwner = -1;
        tetFacei = -1;
        tetPti = -1;

        // Dummy position
        position = pTraits<vector>::max;
    }

    return facei;
}


Foam::label Foam::patchInjectionBaseRPF::setPositionAndCell
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


Foam::label Foam::patchInjectionBaseRPF::whichProc(const scalar fraction01) const
{
    const scalar areaFraction = fraction01*patchArea_;

    // Determine which processor to inject from
    forAllReverse(sumTriMagSf_, i)
    {
        if (areaFraction >= sumTriMagSf_[i])
        {
            return i;


          //  Info <<" processor i    "<< i << endl;
        }
    }

    return 0;
}


// ************************************************************************* //
