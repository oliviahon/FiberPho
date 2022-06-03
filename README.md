# FiberPho
For analysis of fiber photometry recordings and alignment with SimBA output in MATLAB

Written for photometry data recorded with Neurophotometrics FP3001 system. Experiment details here: https://www.biorxiv.org/content/10.1101/2022.03.28.486055v1
additional packages are needed to run initial photometry processing code. Details in FearLearning3ROIregression_GCaMP. 
Behavior videos were tracked post-hoc using DeepLabCut https://github.com/DeepLabCut/DeepLabCut
Freezing behavior was auto-scored using SimBA https://github.com/sgoldenlab/simba

Analysis pipeline is as follows:
1. photometry data initially processed using **FearLearning3ROIregression_GCaMP**, which writes out an excel file with processed, z-scored data for each mouse. 
2. Spikes were detected in cell body recordings using **SpikeThresholding_GCaMP**.
3. Start and end time of each freezing bout was determined using **DetectBouts**
4. Spikes and behavior data were aligned using **SpikesDuringFreezing_GCaMP**. For GCaMP terminal recordings, the average z-score during each bout of freezing and mobility was calculated using **AVGZDuringFreezingEpoch_GCaMP**
