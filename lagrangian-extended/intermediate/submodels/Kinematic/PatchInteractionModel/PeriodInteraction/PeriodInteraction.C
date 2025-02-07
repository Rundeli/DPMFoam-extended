/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2020-2022 OpenCFD Ltd.
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

#include "PeriodInteraction.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class CloudType>
void Foam::PeriodInteraction<CloudType>::writeFileHeader(Ostream& os)
{
    PatchInteractionModel<CloudType>::writeFileHeader(os);

    forAll(nRemovedOnFirstPatch_, i)
    {
        const word& firstPatchName = recyclePatches_[i].first();

        forAll(nRemovedOnFirstPatch_[i], injectori)
        {
            const word suffix = Foam::name(injectori);
            this->writeTabbed(os, firstPatchName + "_nRemoved_" + suffix);
            this->writeTabbed(os, firstPatchName + "_massRemoved_" + suffix);
            this->writeTabbed(os, firstPatchName + "_nInjected_" + suffix);
            this->writeTabbed(os, firstPatchName + "_massInjected_" + suffix);
        }

        const word& secondPatchName = recyclePatches_[i].second();

        forAll(nInjectedOnSecondPatch_[i], injectori)
        {
            const word suffix = Foam::name(injectori);
            this->writeTabbed(os, secondPatchName + "_nRemoved_" + suffix);
            this->writeTabbed(os, secondPatchName + "_massRemoved_" + suffix);
            this->writeTabbed(os, secondPatchName + "_nInjected_" + suffix);
            this->writeTabbed(os, secondPatchName + "_massInjected_" + suffix);
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::PeriodInteraction<CloudType>::PeriodInteraction
(
    const dictionary& dict,
    CloudType& cloud
)
:
    PatchInteractionModel<CloudType>(dict, cloud, typeName),
    mesh_(cloud.mesh()),
    recyclePatches_(this->coeffDict().lookup("recyclePatches")),
    //recyclePatches_.size()返回有几个回收patch pair
    recyclePatchesIds_(recyclePatches_.size()),
    recycledParcels1_(recyclePatches_.size()),
    recycledParcels2_(recyclePatches_.size()),
    //有几对回收patch pair
    nRemovedOnFirstPatch_(recyclePatches_.size()), // per patch the no. of parcels (第一个patch)
    massRemovedOnFirstPatch_(recyclePatches_.size()),
    nInjectedOnFirstPatch_(recyclePatches_.size()),
    massInjectedOnFirstPatch_(recyclePatches_.size()),

    nRemovedOnSecondPatch_(recyclePatches_.size()), //(第二个patch)
    massRemovedOnSecondPatch_(recyclePatches_.size()),
    nInjectedOnSecondPatch_(recyclePatches_.size()),
    massInjectedOnSecondPatch_(recyclePatches_.size()),

    injectionPatchPtr1_(recyclePatches_.size()),//每一对回收patch pair应有两个粒子注入模型
    injectionPatchPtr2_(recyclePatches_.size()),
    recycleFraction_
    (
        this->coeffDict().template getCheck<scalar>
        (
            "recycleFraction",
            scalarMinMax::zero_one()
        )
    ),
    outputByInjectorId_
    (
        this->coeffDict().getOrDefault("outputByInjectorId", false)
    )
{
    // Determine the number of injectors and the injector mapping
    label nInjectors = 0;
    if (outputByInjectorId_)
    {
        for (const auto& inj : cloud.injectors())
        {
            injIdToIndex_.insert(inj.injectorID(), ++nInjectors);
        }
    }

    // The normal case, and safety if injector mapping was somehow null
    if (injIdToIndex_.empty())
    {
        nInjectors = 1;
    }

    forAll(recyclePatches_, i)
    {
        injectionPatchPtr1_.set
        (
            i,
            new PeriodPatchInjectionBase(mesh_, recyclePatches_[i].first(),true,false)
        );
        
        injectionPatchPtr2_.set
        (
            i,
            new PeriodPatchInjectionBase(mesh_, recyclePatches_[i].second(),false,true)
        );

        // Mappings from patch names to patch IDs
        recyclePatchesIds_[i].first() =
            mesh_.boundaryMesh().findPatchID(recyclePatches_[i].first());
        recyclePatchesIds_[i].second() =
            mesh_.boundaryMesh().findPatchID(recyclePatches_[i].second());

        // Storage for reporting
        nRemovedOnFirstPatch_[i].setSize(nInjectors, Zero);
        massRemovedOnFirstPatch_[i].setSize(nInjectors, Zero);
        nInjectedOnFirstPatch_[i].setSize(nInjectors, Zero);
        massInjectedOnFirstPatch_[i].setSize(nInjectors, Zero);

        nRemovedOnSecondPatch_[i].setSize(nInjectors, Zero);
        massRemovedOnSecondPatch_[i].setSize(nInjectors, Zero);
        nInjectedOnSecondPatch_[i].setSize(nInjectors, Zero);
        massInjectedOnSecondPatch_[i].setSize(nInjectors, Zero);
    }

}


template<class CloudType>
Foam::PeriodInteraction<CloudType>::PeriodInteraction
(
    const PeriodInteraction<CloudType>& pim
)
:
    PatchInteractionModel<CloudType>(pim),
    mesh_(pim.mesh_),
    recyclePatches_(pim.recyclePatches_),
    recyclePatchesIds_(pim.recyclePatchesIds_),
    recycledParcels1_(pim.recycledParcels1_),
    recycledParcels2_(pim.recycledParcels2_),
    nRemovedOnFirstPatch_(pim.nRemovedOnFirstPatch_),
    massRemovedOnFirstPatch_(pim.massRemovedOnFirstPatch_),
    nInjectedOnFirstPatch_(pim.nInjectedOnFirstPatch_),
    massInjectedOnFirstPatch_(pim.massInjectedOnFirstPatch_),
    nRemovedOnSecondPatch_(pim.nRemovedOnSecondPatch_), 
    massRemovedOnSecondPatch_(pim.massRemovedOnSecondPatch_),
    nInjectedOnSecondPatch_(pim.nInjectedOnSecondPatch_),
    massInjectedOnSecondPatch_(pim.massInjectedOnSecondPatch_),
    injIdToIndex_(pim.injIdToIndex_),
    recycleFraction_(pim.recycleFraction_),
    outputByInjectorId_(pim.outputByInjectorId_)
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::PeriodInteraction<CloudType>::correct
(
    typename CloudType::parcelType& p,
    const polyPatch& pp,
    bool& keepParticle
)
{
    // Injector ID
    const label idx =
    (
        injIdToIndex_.size()
      ? injIdToIndex_.lookup(p.typeId(), 0)
      : 0
    );

    //检查该patch是否在recyclePatch中
    label addri = -1;
    bool isFirstPatch = false;
    bool isSecondPatch = false;
    forAll(recyclePatchesIds_, i)
    {
        if (recyclePatchesIds_[i].contains(pp.index())) //该patch位于回收列表中
        {
            addri = i;

            if (recyclePatchesIds_[i].first() == pp.index()) //该patch是回收pair中第一个patch
            {
                isFirstPatch = true;
                break;
            }
            else //该patch是回收pair中第二个patch
            {
                isSecondPatch = true;
                break;
            }
        }
        else //该patch不在回收patch pair中
        {
            keepParticle = true;
            return false;
        }
    }

    // Flag to remove current parcel and copy to local storage
    keepParticle = false;

    if(isFirstPatch) //第一个patch上的粒子回收至第一个patch上的recycledParcels1_中
    {
        recycledParcels1_[addri].append
        (
            static_cast<parcelType*>(p.clone().ptr())
        );

        ++nRemovedOnFirstPatch_[addri][idx];
        massRemovedOnFirstPatch_[addri][idx] += p.nParticle()*p.mass();
    }
    if(isSecondPatch) //第二个patch上的粒子回收至第二个patch上的recycledParcels2_中
    {
        recycledParcels2_[addri].append
        (
            static_cast<parcelType*>(p.clone().ptr())
        );

        ++nRemovedOnSecondPatch_[addri][idx];
        massRemovedOnSecondPatch_[addri][idx] += p.nParticle()*p.mass();
    }

    return true;
}


template<class CloudType>
void Foam::PeriodInteraction<CloudType>::postEvolve()
{

    if (Pstream::parRun()) //并行运行
    {
        // See comments in Cloud::move() about transfer particles handling

        // Allocate transfer buffers
        PstreamBuffers pBufs1;
        PstreamBuffers pBufs2;

        // Cache of opened UOPstream wrappers
        PtrList<UOPstream> UOPstreamPtrs(Pstream::nProcs());

        forAll(recycledParcels1_, addri) //对所有回收patch pair中第一个patch所需回收粒子进行遍历
        {
            auto& patchParcels = recycledParcels1_[addri];
            //auto& injectionPatch = injectionPatchPtr_[addri];

            for (parcelType& p : patchParcels) //对第一个patch上的待回收粒子进行处理
            {
                point positionOnExit = p.position();

                // Identify the processor that owns the location
                //显示为硬编码为0，目前不能对patch进行分解分属多个processor，下一步研究如何计算processor编号方法（该作用域无法访问Triface_计算粒子注入位置）
                const label toProci = 5;

                // Get/create output stream
                auto* osptr = UOPstreamPtrs.get(toProci);
                if (!osptr)
                {
                    osptr = new UOPstream(toProci, pBufs1);
                    UOPstreamPtrs.set(toProci, osptr);
                }

                // Tuple: (address positionOnExit particle)
                (*osptr) << addri << positionOnExit << p;

                // Can now remove from list and delete
                delete(patchParcels.remove(&p));

            }
        }

        forAll(recycledParcels2_, addri) //对所有回收patch pair中第二个patch所需回收粒子进行遍历
        {
            auto& patchParcels = recycledParcels2_[addri];
            
            for (parcelType& p : patchParcels) //对第二个patch上的待回收粒子进行处理
            {
                point positionOnExit = p.position();

                // Identify the processor that owns the location
                //显示为硬编码0，目前不能对patch进行分解分属多个processor，下一步研究如何计算processor编号方法（该作用域无法访问Triface_计算粒子注入位置）
                const label toProci = 0;

                // Get/create output stream
                auto* osptr = UOPstreamPtrs.get(toProci);
                if (!osptr)
                {
                    osptr = new UOPstream(toProci, pBufs2);
                    UOPstreamPtrs.set(toProci, osptr);
                }

                // Tuple: (address positionOnExit particle)
                (*osptr) << addri << positionOnExit << p;

                // Can now remove from list and delete
                delete(patchParcels.remove(&p));
            }
        }

        pBufs1.finishedSends();
        pBufs2.finishedSends();

        // Not looping, so no early exit needed
        //
        // if (!returnReduceOr(pBufs.hasRecvData()))
        // {
        //     // No parcels to recycle
        //     return;
        // }

        // Retrieve from receive buffers
        for (const int proci : pBufs1.allProcs())
        {
            if (pBufs1.recvDataCount(proci))
            {

                UIPstream is(proci, pBufs1);
                
                // Read out each (address positionOnExit particle) tuple
                while (!is.eof())
                {
                    const label addri = pTraits<label>(is);
                    const point positionOnExit = pTraits<vector>(is);
                    auto* newp = new parcelType(this->owner().mesh(), is);

                    point positionOnEntry;
                    label enterTrifaceIndex;
                    label cellOwner;
                    label dummy;

                    //默认第一个面处于processor0，将第二个patch上的回收粒子位置数据传输至processor0，调用第一个面上的粒子注入模型injectionPatchPtr1_[addri]
                    
                    const vector direction(1,0,0); //第二个patch朝向的第一个patch的垂向方向
                    injectionPatchPtr2_[addri].setEntryPos(mesh_,positionOnExit,direction,positionOnEntry,enterTrifaceIndex);

                    // Parcel to be recycled
                    injectionPatchPtr2_[addri].setPositionAndCell
                    (
                        mesh_,
                        this->owner().rndGen(),
                        positionOnEntry,
                        enterTrifaceIndex,
                        cellOwner,
                        dummy,
                        dummy
                    );

                    const label idx =
                    (
                        injIdToIndex_.size()
                        ? injIdToIndex_.lookup(newp->typeId(), 0)
                        : 0
                    );

                    //统计数据
                    ++nInjectedOnSecondPatch_[addri][idx];
                    massInjectedOnSecondPatch_[addri][idx] += newp->nParticle()*newp->mass();

                    //注入粒子
                    newp->relocate(positionOnEntry,cellOwner);
                    this->owner().addParticle(newp);
                    
                }
            }
        }

        for (const int proci : pBufs2.allProcs())
        {
            if (pBufs2.recvDataCount(proci))
            {
                UIPstream is(proci, pBufs2);
                
                // Read out each (address positionOnExit particle) tuple
                while (!is.eof())
                {
                    const label addri = pTraits<label>(is);
                    const point positionOnExit = pTraits<vector>(is);
                    auto* newp = new parcelType(this->owner().mesh(), is);

                    point positionOnEntry;
                    label enterTrifaceIndex;
                    label cellOwner;
                    label dummy;

                    const vector direction(-1,0,0); //第一个patch朝向第二个patch的垂向方向
                    injectionPatchPtr1_[addri].setEntryPos(mesh_,positionOnExit,direction,positionOnEntry,enterTrifaceIndex);
                        
                    // Parcel to be recycled
                    injectionPatchPtr1_[addri].setPositionAndCell
                    (
                        mesh_,
                        this->owner().rndGen(),
                        positionOnEntry,
                        enterTrifaceIndex,
                        cellOwner,
                        dummy,
                        dummy
                    );

                    const label idx =
                    (
                        injIdToIndex_.size()
                        ? injIdToIndex_.lookup(newp->typeId(), 0)
                        : 0
                    );

                    //统计数据
                    ++nInjectedOnFirstPatch_[addri][idx];
                    massInjectedOnFirstPatch_[addri][idx] += newp->nParticle()*newp->mass();

                    //注入粒子
                    newp->relocate(positionOnEntry,cellOwner);
                    this->owner().addParticle(newp);
                }
            }
        }
    }
    else //单核运行
    {
        forAll(recycledParcels1_, addri) //对于每对回收patch pair中第一个patch回收粒子列表进行遍历
        {
            for (parcelType& p : recycledParcels1_[addri]) //遍历其中一对回收patch pair中第一个patch回收粒子列表
            {
                //获取旧粒子位置坐标
                point positionOnExit = p.position();
                
                //粒子移出回收粒子list 
                parcelType* newp = recycledParcels1_[addri].remove(&p);

                //初始化其余setPositionAndCell()需传入参数
                point positionOnEntry;
                label enterTrifaceIndex;
                label cellOwner;
                label dummy;
                const vector direction(1,0,0);

                //第一个patch上的回收粒子应在第二个patch中注入，使用第二个patch的粒子注入模型injectionPatchPtr2_
                injectionPatchPtr2_[addri].setPositionAndCell
                (
                    mesh_,
                    this->owner().rndGen(),
                    positionOnExit,
                    positionOnEntry,
                    direction,
                    enterTrifaceIndex,
                    cellOwner,
                    dummy,
                    dummy
                );
        
                //调用particle::relocate()移动粒子
                newp->relocate(positionOnEntry,cellOwner);

                // Injector ID
                const label idx =
                (
                    injIdToIndex_.size()
                  ? injIdToIndex_.lookup(newp->typeId(), 0)
                  : 0
                );
                ++nInjectedOnSecondPatch_[addri][idx];
                massInjectedOnSecondPatch_[addri][idx] += newp->nParticle()*newp->mass();

                this->owner().addParticle(newp);
            }
        }

        forAll(recycledParcels2_, addri)// 对于每对回收patch pair中第一个patch回收粒子列表进行遍历
        {
            for (parcelType& p : recycledParcels2_[addri])
            {
                //获取旧粒子位置坐标
                point positionOnExit = p.position();
                
                //粒子移出回收粒子list
                parcelType* newp = recycledParcels2_[addri].remove(&p);

                //初始化其余setPositionAndCell()需传入参数
                point positionOnEntry;
                label enterTrifaceIndex;
                label cellOwner;
                label dummy;
                const vector direction(-1,0,0);

                //调用PeriodPatchInjectionBase::setPositionAndCell()中
                injectionPatchPtr1_[addri].setPositionAndCell
                (
                    mesh_,
                    this->owner().rndGen(),
                    positionOnExit,
                    positionOnEntry,
                    direction,
                    enterTrifaceIndex,
                    cellOwner,
                    dummy,
                    dummy
                );
        
                //调用particle::relocate()移动粒子
                newp->relocate(positionOnEntry,cellOwner);

                // Injector ID
                const label idx =
                (
                    injIdToIndex_.size()
                  ? injIdToIndex_.lookup(newp->typeId(), 0)
                  : 0
                );
                ++nInjectedOnFirstPatch_[addri][idx];
                massInjectedOnFirstPatch_[addri][idx] += newp->nParticle()*newp->mass();

                this->owner().addParticle(newp);
            }
        }
    }
}


template<class CloudType>
void Foam::PeriodInteraction<CloudType>::info()
{
    //提前声明包含PatchInteractionModel中的info函数
    PatchInteractionModel<CloudType>::info();


    labelListList npr0(nRemovedOnFirstPatch_.size());
    scalarListList mpr0(massRemovedOnFirstPatch_.size());
    labelListList npi0(nInjectedOnFirstPatch_.size());
    scalarListList mpi0(massInjectedOnFirstPatch_.size());

    labelListList npr1(nRemovedOnSecondPatch_.size());
    scalarListList mpr1(massRemovedOnSecondPatch_.size());
    labelListList npi1(nInjectedOnSecondPatch_.size());
    scalarListList mpi1(massInjectedOnSecondPatch_.size());

    forAll(nRemovedOnFirstPatch_, patchi)
    {
        const label lsd = nRemovedOnFirstPatch_[patchi].size();
        npr0[patchi].setSize(lsd, Zero);
        mpr0[patchi].setSize(lsd, Zero);
        npi0[patchi].setSize(lsd, Zero);
        mpi0[patchi].setSize(lsd, Zero);
    }

    forAll(nRemovedOnSecondPatch_, patchi)
    {
        const label lsd = nRemovedOnSecondPatch_[patchi].size();
        npr1[patchi].setSize(lsd, Zero);
        mpr1[patchi].setSize(lsd, Zero);
        npi1[patchi].setSize(lsd, Zero);
        mpi1[patchi].setSize(lsd, Zero);
    }

    this->getModelProperty("nRemovedOnFirstPatch", npr0);
    this->getModelProperty("massRemovedOnFirstPatch", mpr0);
    this->getModelProperty("nInjectedOnFirstPatch", npi0);
    this->getModelProperty("massInjectedOnFirstPatch", mpi0);

    this->getModelProperty("nRemovedOnSecondPatch", npr1);
    this->getModelProperty("massRemovedOnSecondPatch", mpr1);
    this->getModelProperty("nInjectedOnSecondPatch", npi1);
    this->getModelProperty("massInjectedOnSecondPatch", mpi1);

    // Accumulate current data
    labelListList nprf(nRemovedOnFirstPatch_);

    forAll(nprf, i)         //for first patch
    {
        Pstream::listCombineGather(nprf[i], plusEqOp<label>());
        nprf[i] = nprf[i] + npr0[i];
    }

    scalarListList mprf(massRemovedOnFirstPatch_);
    forAll(mprf, i)
    {
        Pstream::listCombineGather(mprf[i], plusEqOp<scalar>());
        mprf[i] = mprf[i] + mpr0[i];
    }

    labelListList npif(nInjectedOnFirstPatch_);
    forAll(npif, i)
    {
        Pstream::listCombineGather(npif[i], plusEqOp<label>());
        npif[i] = npif[i] + npi0[i];
    }

    scalarListList mpif(massInjectedOnFirstPatch_);
    forAll(mpif, i)
    {
        Pstream::listCombineGather(mpif[i], plusEqOp<scalar>());
        mpif[i] = mpif[i] + mpi0[i];
    }

    labelListList nprs(nRemovedOnSecondPatch_);
    forAll(nprs, i)     //for second patch
    {
        Pstream::listCombineGather(nprs[i], plusEqOp<label>());
        nprs[i] = nprs[i] + npr1[i];
    }

    scalarListList mprs(massRemovedOnSecondPatch_);
    forAll(mprs, i)
    {
        Pstream::listCombineGather(mprs[i], plusEqOp<scalar>());
        mprs[i] = mprs[i] + mpr1[i];
    }

    labelListList npis(nInjectedOnSecondPatch_);
    forAll(npis, i)
    {
        Pstream::listCombineGather(npis[i], plusEqOp<label>());
        npis[i] = npis[i] + npi1[i];
    }

    scalarListList mpis(massInjectedOnSecondPatch_);
    forAll(mpis, i)
    {
        Pstream::listCombineGather(mpis[i], plusEqOp<scalar>());
        mpis[i] = mpis[i] + mpi1[i];
    }

    if (injIdToIndex_.size())
    {
        // Since injIdToIndex_ is a one-to-one mapping (starting as zero),
        // can simply invert it.
        labelList indexToInjector(injIdToIndex_.size());
        forAllConstIters(injIdToIndex_, iter)
        {
            indexToInjector[iter.val()] = iter.key();
        }

        forAll(nprf, i)
        {
            const word& firstPatchName =  recyclePatches_[i].first();

            Log_<< "    Parcel fate: patch " <<  firstPatchName
                << " (number, mass)" << nl;

            forAll(mprf[i], indexi)
            {
                Log_<< "      - removed  (injector " << indexToInjector[indexi]
                    << ")  = " << nprf[i][indexi]
                    << ", " << mprf[i][indexi] << nl
                    << "      - injected  (injector " << indexToInjector[indexi]
                    << ")  = " << npif[i][indexi]
                    << ", " << mpif[i][indexi] << nl;
                this->file()
                    << tab << nprf[i][indexi] << tab << mprf[i][indexi]
                    << tab << npif[i][indexi] << tab << mpif[i][indexi];
            }

            const word& secondPatchName = recyclePatches_[i].second();

            Log_<< "    Parcel fate: patch " <<  secondPatchName
                << " (number, mass)" << nl;

            forAll(mpis[i], indexi)
            {
                Log_<< "      - removed  (injector " << indexToInjector[indexi]
                    << ")  = " << nprs[i][indexi]
                    << ", " << mprs[i][indexi] << nl
                    << "      - injected  (injector " << indexToInjector[indexi]
                    << ")  = " << npis[i][indexi]
                    << ", " << mpis[i][indexi] << nl;
                this->file()
                    << tab << nprs[i][indexi] << tab << mprs[i][indexi]
                    << tab << npis[i][indexi] << tab << mpis[i][indexi];
            }
        }

        this->file() << endl;
    }
    else
    {
        forAll(nprf, i)
        {
            const word& firstPatchName =  recyclePatches_[i].first();

            Log_<< "    Parcel fate: patch " <<  firstPatchName
                << " (number, mass)" << nl
                << "      - removed    = " << nprf[i][0] << ", " << mprf[i][0]
                << nl;

            Log_<< "      - injected    = " << npif[i][0] << ", " << mpif[i][0]
                << nl;

            this->file()
                << tab << nprf[i][0] << tab << mprf[i][0] << tab << npif[i][0] << tab << mpif[i][0];
        }

        forAll(npis, i)
        {
            const word& secondPatchName = recyclePatches_[i].second();

            Log_<< "    Parcel fate: patch " <<  secondPatchName
                << " (number, mass)" << nl
                << "      - removed    = " << nprs[i][0] << ", " << mprs[i][0]
                << nl;

            Log_<< "      - injected   = " << npis[i][0] << ", " << mpis[i][0]
                << nl;

            this->file()
                << tab << nprs[i][0] << tab << mprs[i][0] << tab << npis[i][0] << tab << mpis[i][0];
        }

        this->file() << endl;
    }

    if (this->writeTime())
    {
        this->setModelProperty("nRemovedOnFirstPatch", nprf);
        this->setModelProperty("massRemovedOnFirstPatch", mprf);
        this->setModelProperty("nInjectedOnFirstPatch", npif);
        this->setModelProperty("massInjectedOnFirstPatch", mpif);

        this->setModelProperty("nRemovedOnSecondPatch", nprs);
        this->setModelProperty("massRemovedOnSecondPatch", mprs);
        this->setModelProperty("nInjectedOnSecondPatch", npis);
        this->setModelProperty("massInjectedOnSecondPatch", mpis);

        nRemovedOnFirstPatch_ = Zero;
        massRemovedOnFirstPatch_ = Zero;
        nInjectedOnFirstPatch_ = Zero;
        massInjectedOnFirstPatch_ = Zero;

        nRemovedOnSecondPatch_ = Zero;
        massRemovedOnSecondPatch_ = Zero;
        nInjectedOnSecondPatch_ = Zero;
        massInjectedOnSecondPatch_ = Zero;
    }
}


// ************************************************************************* //