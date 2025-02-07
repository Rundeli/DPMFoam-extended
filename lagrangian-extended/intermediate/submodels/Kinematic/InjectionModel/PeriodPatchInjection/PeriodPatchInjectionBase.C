/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2024 OpenCFD Ltd.
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

#include "PeriodPatchInjectionBase.H"
#include "polyMesh.H"
#include "SubField.H"
#include "Random.H"
#include "triangle.H"
#include "volFields.H"
#include "polyMeshTetDecomposition.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::PeriodPatchInjectionBase::PeriodPatchInjectionBase
(
    const polyMesh& mesh,
    const word& patchName,
    const bool& isFirstPatch,
    const bool& isSecondPatch
)
:
    patchName_(patchName),
    patchId_(mesh.boundaryMesh().findPatchID(patchName_)),
    patchArea_(0),
    patchNormal_(),
    cellOwners_(),
    triFace_(),
    triCumulativeMagSf_(),
    sumTriMagSf_(),
    isFirstPatch_(isFirstPatch),
    isSecondPatch_(isSecondPatch)
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


Foam::PeriodPatchInjectionBase::PeriodPatchInjectionBase(const PeriodPatchInjectionBase& pib)
:
    patchName_(pib.patchName_),
    patchId_(pib.patchId_),
    patchArea_(pib.patchArea_),
    patchNormal_(pib.patchNormal_),
    cellOwners_(pib.cellOwners_),
    triFace_(pib.triFace_),
    triCumulativeMagSf_(pib.triCumulativeMagSf_),
    sumTriMagSf_(pib.sumTriMagSf_),
    isFirstPatch_(pib.isFirstPatch_),
    isSecondPatch_(pib.isSecondPatch_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::PeriodPatchInjectionBase::updateMesh(const polyMesh& mesh)
{
    // Set/cache the injector cells
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];
    const pointField& points = patch.points();

    //得到injector patch上所有face所属cells的labelList
    cellOwners_ = patch.faceCells();

    // Triangulate the patch faces and create addressing
    {
        //初始化三角面个数为0
        label nTris = 0;
        //遍历patch上的所有face，计算每个face上分割得到的三角面个数，累加得到patch上总的三角面个数
        for (const face& f : patch)
        {
            nTris += f.nTriangles();
        }

        //初始化一个元素为labelledTri的动态列表，初始长度为patch上三角面个数
        DynamicList<labelledTri> dynTriFace(nTris);

        //初始化一个元素为face的动态列表，初始长度为8
        DynamicList<face> tris(8);  // work array

        //遍历patch上每一个面
        forAll(patch, facei)
        {
            const face& f = patch[facei];

            //清空tris
            tris.clear();
            //将face分割为三角面，将所得三角面存入tris中
            f.triangles(points, tris);

            //遍历tris中所有三角面，以三角面的三个点以及所属face创建labelledTri对象存入dynTriFace中
            for (const auto& t : tris)
            {
                dynTriFace.emplace_back(t[0], t[1], t[2], facei);
            }
        }

        // Transfer to persistent storage
        triFace_.transfer(dynTriFace);

    }

    //初始化当前processor序号
    const label myProci = UPstream::myProcNo();
    //初始化processor总数
    const label numProc = UPstream::nProcs();

    // Calculate the cumulative triangle weights
    {
        //初始化一个长度为TriFace+1的scalarList
        triCumulativeMagSf_.resize_nocopy(triFace_.size()+1);

        //初始化指向scalarList起始位置的迭代器
        auto iter = triCumulativeMagSf_.begin();

        // Set zero value at the start of the tri area/weight list
        scalar patchArea = 0;
        //迭代器向后移动一个元素，将其值设为patchArea
        *iter++ = patchArea;

        // Calculate cumulative and total area
        for (const auto& t : triFace_)
        {   
            //前一个triFace处累积patchArea加上当前TriFace面积
            patchArea += t.mag(points);
            //将当前Triface处triCumulativeMagSf_元素设为累加当前Triface面积后的patchArea
            *iter++ = patchArea;
        }

        //初始化一个长度为processor+1的scalarList
        sumTriMagSf_.resize_nocopy(numProc+1);
        //将第一个元素设为0
        sumTriMagSf_[0] = 0;

        {
            scalarList::subList slice(sumTriMagSf_, numProc, 1);
            slice[myProci] = patchArea;
            Pstream::allGatherList(slice);
        }

        // Convert to cumulative
        for (label i = 1; i < sumTriMagSf_.size(); ++i)
        {
            sumTriMagSf_[i] += sumTriMagSf_[i-1];
        }
    }

    //初始化patch上face面法向向量模的标量场magSf
    const scalarField magSf(mag(patch.faceAreas()));

    //累加得到patch总面积patchArea_
    patchArea_ = sum(magSf);

    //得到face的法向单位向量
    patchNormal_ = patch.faceAreas()/magSf;
    reduce(patchArea_, sumOp<scalar>());
}

//called directly in parallel case with calculated faction01
Foam::label Foam::PeriodPatchInjectionBase::setPositionAndCell
(
    const fvMesh& mesh,
    Random& rnd,
    vector& position,
    const label& trii,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{

    //初始化face索引为-1
    label facei = -1;

    //如果face所属cell列表不为空
    if (!cellOwners_.empty())
    {

        label proci = -1;

        // Determine which processor to inject from
        if(isFirstPatch_)
        {
            proci = 0;
        }
        else if(isSecondPatch_)
        {
            proci = 5;
        }
        else
        {
            FatalErrorInFunction
                << "Must inject from the the first patch or the second patch, check your initalization statement"
                << exit(FatalError);
        }

        //如果当前processor是粒子注入位置所属processor
        if (UPstream::myProcNo() == proci)
        {
            // Set cellOwner
            facei = triFace_[trii].index();
            cellOwner = cellOwners_[trii];

            vector positionBk = position;

            // Position perturbed away from face (into domain)
            const scalar a = rnd.position(scalar(0.1), scalar(0.5));
            const vector& pc = mesh.cellCentres()[cellOwner];
            const vector d =
                mag((positionBk - pc) & patchNormal_[facei])*patchNormal_[facei];

            position = positionBk - 0.01*a*d;

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
    }

    if (facei == -1)
    {
        cellOwner = -1;
        tetFacei = -1;
        tetPti = -1;

        // Dummy position
        position = pTraits<vector>::max;
    }

    return facei;
}

//called in single-core case calculate fraction01 after function being called
Foam::label Foam::PeriodPatchInjectionBase::setPositionAndCell
(
    const fvMesh& mesh,
    Random& rnd,
    point& position,
    point& positionOnEntry,
    const vector& direction,
    label& trii,
    label& cellOwner,
    label& tetFacei,
    label& tetPti
)
{
    label enterTrifaceIndex;
    setEntryPos(mesh,position,direction,positionOnEntry,enterTrifaceIndex);

    return setPositionAndCell
    (
        mesh,
        rnd,
        positionOnEntry,
        enterTrifaceIndex,
        cellOwner,
        tetFacei,
        tetPti
    );
}
 
void Foam::PeriodPatchInjectionBase::setEntryPos
( 
    const fvMesh& mesh,
    const point& positionOnExit,
    const vector& direction,
    vector& positionOnEntry,
    label& enterTrifaceIndex
)
{
    //初始化边界面patch
    const polyPatch& patch = mesh.boundaryMesh()[patchId_];
    //初始化该patch点集
    const pointField& points = patch.points();

    scalar i = 0; 
    enterTrifaceIndex = -1;
    
    for (const auto& t : triFace_)
    {
        pointHit pHit = t.intersection(positionOnExit,direction,points,intersection::HALF_RAY);
        if(pHit.hit())
        {
            positionOnEntry = pHit.hitPoint();
            enterTrifaceIndex = i;
            break;
        }
        i++;
    }

    if (enterTrifaceIndex == -1)
    {
        FatalErrorInFunction
            << "The search for entry location failed. On patch: "
            << patchName_ 
            << positionOnExit
            << positionOnEntry
            << enterTrifaceIndex <<" "
            << i
            << exit(FatalError);
    }
}


/*Foam::label Foam::PeriodPatchInjectionBase::whichProc(const scalar fraction01) const
{


    const scalar areaFraction = fraction01*patchArea_;

    // Determine which processor to inject from
    forAllReverse(sumTriMagSf_, i)
    {
        if (areaFraction >= sumTriMagSf_[i])
        {
            return i;
        }
    }

    return 0;
}*/
// ************************************************************************* //