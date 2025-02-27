/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::PeriodInteraction

Group
    grpLagrangianIntermediatePatchInteractionSubModels

Description
    Patch interaction model to perform
    parcel transfer and recycle from one patch to another.

    

Usage
    Minimal example by using
    \c constant/reactingCloud1Properties.subModels.multiInteractionCoeffs
    \verbatim
    multiInteractionCoeffs
    {
        ...

        model1
        {
            // Mandatory entries (unmodifiable)
            patchInteractionModel    periodInteraction;

            periodInteractionCoeffs
            {
                recyclePatches
                (
                    (<outletPatch1> <inletPatch1>)
                    (<outletPatch2> <inletPatch2>)
                    ...
                );

                recycleFraction      0.8; // [0-1]

                // Optional entries (unmodifiable)
                outputByInjectorId   false;
            }
        }

        ...
    }
    \endverbatim

    where the entries mean:
    \table
      Property     | Description                      | Type | Reqd | Dflt
      patchInteractionModel | Type name: PeriodInteraction  <!--
                                                  --> | word | yes  | -
      recyclePatches | Names of outlet-inlet patch pairs     <!--
                                                -->   | (word word) | yes | -
      recycleFraction | Fraction of parcels recycled from <!--
                                       --> outlet to inlet | scalar | yes | -
      outputByInjectorId | Flag to output escaped/mass    <!--
                     --> particles sorted by injectorID | bool | no | false
    \endtable

SourceFiles
    PeriodInteraction.C

\*---------------------------------------------------------------------------*/

#ifndef PeriodInteraction_H
#define PeriodInteraction_H

#include "PeriodPatchInjectionBase.H"
#include "IDLList.H"
#include "PtrList.H"
#include "PatchInteractionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class PeriodInteraction Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class PeriodInteraction
:
    public PatchInteractionModel<CloudType>
{
protected:

    typedef typename CloudType::parcelType parcelType;

    // Protected Data

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Outlet-inlet patch pair to apply parcel recycling
        List<Pair<word>> recyclePatches_;

        //- Patch IDs of recyclePatches
        List<Pair<label>> recyclePatchesIds_;

        //- Parcel IDs of recycled parcels
        List<IDLList<parcelType>> recycledParcels1_;

        List<IDLList<parcelType>> recycledParcels2_;


        // Bookkeeping for particle fates

            //- Number of parcels removed
            List<List<label>> nRemovedOnFirstPatch_;

            //- Mass of parcels removed
            List<List<scalar>> massRemovedOnFirstPatch_;

            //- Number of parcels injected
            List<List<label>> nInjectedOnFirstPatch_;

            //- Mass of parcels injected
            List<List<scalar>> massInjectedOnFirstPatch_;

            //第二个patch上的统计数据
            List<List<label>> nRemovedOnSecondPatch_;

            List<List<scalar>> massRemovedOnSecondPatch_;

            List<List<label>> nInjectedOnSecondPatch_;

            List<List<scalar>> massInjectedOnSecondPatch_;

        //- Injector ID to local index map
        Map<label> injIdToIndex_;

        //- Injection patch pointer //双向周期性回收模型
        PtrList<PeriodPatchInjectionBase> injectionPatchPtr1_;

        PtrList<PeriodPatchInjectionBase> injectionPatchPtr2_;
        // List<PtrList<PeriodPatchInjectionBase>> injectionPatchPtr_;

        //- Parcel fraction to be recycled from outlet to inlet
        const scalar recycleFraction_;

        //- Flag to output escaped/mass particles sorted by injectorID
        bool outputByInjectorId_;

    // Protected Member Functions

        //- Output file header information
        virtual void writeFileHeader(Ostream& os);


public:

    //- Runtime type information
    TypeName("periodInteraction");

    // Constructors

        //- Construct from dictionary
        PeriodInteraction(const dictionary& dict, CloudType& cloud);

        //- Construct copy from owner cloud and patch interaction model
        PeriodInteraction(const PeriodInteraction<CloudType>& pim);

        //- Construct and return a clone using supplied owner cloud
        virtual autoPtr<PatchInteractionModel<CloudType>> clone() const
        {
            return autoPtr<PatchInteractionModel<CloudType>>
            (
                new PeriodInteraction<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~PeriodInteraction() = default;


    // Member Functions

        //- Apply velocity correction
        //  Returns true if particle remains in same cell
        virtual bool correct
        (
            typename CloudType::parcelType& p,
            const polyPatch& pp,
            bool& keepParticle
        );

        //- Post-evolve hook
        virtual void postEvolve();

        //- Write patch interaction info
        virtual void info();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "PeriodInteraction.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //