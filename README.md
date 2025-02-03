这是一个与原OpenFOAM lagrangrian库代码进行隔离，以便自由扩展的DPMFoam版本，目前扩展工作有：
1.新增单向周期性粒子patchInteractionModel: periodInteraction，粒子在outflow patch离开，以离开时在patch上相同位置与速度注入inflow patch，用以模拟无限长的周期性物理条件相似水槽，避免粒子在模拟过程中流失。

下一步工作:
1.改为双向周期性模型，以模拟波浪等往复流作用下无限上的周期性物理条件水槽
2.在DPMFoam中新增对用户自定义源项的支持（体积力）

This is a extended version of DPMFoam which is seperate from orginal OpenFOAM lagrangrian code to aviod compiling problem when implementing new patchInteraction models, particle models and so on.
The following model has been implemented:
1.One-direction particle period recycle patchInteractionModel: periodInteraction.Particles exit the outflow patch and are injected into the inflow patch with the same position and velocity they had upon exiting the patch，which is intened to simulate a wave flume which has infinite length and aviod particle escape during simulation.

Things to do next:
1:add support for two-direction period recycle to simulate wave flume under periodic wave;
2.add support for custom source term in DPMFoam (wave Volumetric force).
