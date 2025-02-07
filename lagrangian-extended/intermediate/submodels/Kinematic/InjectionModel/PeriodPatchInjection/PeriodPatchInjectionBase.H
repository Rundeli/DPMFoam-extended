/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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

Class
    Foam::PatchInjectionBase

Description
    Base class for patch-based injection models.

    Class handles injecting at a random point adjacent to the patch faces to
    provide a more stochastic view of the injection process.  Patch faces are
    triangulated, and area fractions accumulated.  The fractional areas are
    then applied to determine across which face a parcel is to be injected.

SourceFiles
    PatchInjectionBase.C

\*---------------------------------------------------------------------------*/

#ifndef Foam_PeriodPatchInjectionBase_H
#define Foam_PeriodPatchInjectionBase_H

#include "word.H"
#include "labelList.H"
#include "scalarList.H"
#include "vectorList.H"
#include "faceList.H"
#include "labelledTri.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward Declarations
class polyMesh;
class fvMesh;
class Random;

/*---------------------------------------------------------------------------*\
                     Class PeriodPatchInjectionBase Declaration
\*---------------------------------------------------------------------------*/

class PeriodPatchInjectionBase
{
protected:

    // Protected Data

        //- Patch name
        const word patchName_;

        //- Patch ID
        const label patchId_;

        //- Patch area - total across all processors
        scalar patchArea_;

        //- Patch face normal directions
        vectorList patchNormal_;

        //- List of cell labels corresponding to injector positions
        labelList cellOwners_;

        //- The polyPatch faces as triangles, the index of each corresponds
        //- to the undecomposed patch face index.
        List<labelledTri> triFace_;

        //- Cumulative triangle area per triangle face (processor-local)
        scalarList triCumulativeMagSf_;

        //- Cumulative area fractions per processor
        scalarList sumTriMagSf_;

        //定位粒子注入模型所属patch，以确定所属processor
        bool isFirstPatch_;

        bool isSecondPatch_;

public:

    // Constructors

        //- Construct from mesh and patch name
        PeriodPatchInjectionBase(const polyMesh& mesh, const word& patchName, const bool& isFirstPatch, const bool& isSecondPatch);

        //- Copy constructor
        PeriodPatchInjectionBase(const PeriodPatchInjectionBase& pib);


    //- Destructor
    virtual ~PeriodPatchInjectionBase() = default;


    // Member Functions

        //- Update patch geometry and derived info for injection locations
        virtual void updateMesh(const polyMesh& mesh);

        //- Set the injection position and owner cell, tetFace and tetPt
        //  Supply the fraction used to determine the location on the patch
        //  Returns the seed patch face index
        label setPositionAndCell
        (
            const fvMesh& mesh,
            Random& rnd,
            vector& position,
            const label& trii,
            label& cellOwner,
            label& tetFacei,
            label& tetPti
        );

        //- Set the injection position and owner cell, tetFace and tetPt
        //  Returns the seed patch face index
        virtual Foam::label setPositionAndCell
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
        );

        void setEntryPos
        (
            const fvMesh& mesh,
            const point& positionOnExit,
            const vector& direction,
            point& positionOnEntry,
            label& enterTrifaceIndex
        );

        //- Return the processor that has the location specified by the fraction
        label whichProc(const scalar fraction01) const;
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //