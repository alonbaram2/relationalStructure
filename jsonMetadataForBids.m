addpath(genpath('/home/fs0/abaram/scratch/MATLAB/jsonlab-1.5'))

subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07',...
    'sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19',...
    'sub-20','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28'};

opt.ForceRootName = 0; % root name kept empty

% metadata for all BOLD blocks.
objBold.RepetitionTime = 1.235;
objBold.EchoTime = 0.020;
objBold.FlipAngle = 65;
objBold.MultibandAccelerationFactor = 3;
objBold.ParallelReductionFactorInPlane = 2;
objBold.PhaseEncodingDirection = 'j-';
for iSub=1:length(subjects)
    for r=1:2
        for c = 1:4
            opt.FileName = ['/vols/Scratch/abaram/NN2020_BIDS/' subjects{iSub} '/func/' subjects{iSub} '_task-cond0' num2str(c) '_run-0' num2str(r) '_bold.json'];
            objBold.TaskName = ['cond0' int2str(c)];
            savejson('',objBold,opt);
            
        end
    end
end
%%
% fieldmap
objPhase.EchoTime1 = 0.00492;
objPhase.EchoTime2 = 0.00738;
for iSub=1:length(subjects)
    opt.FileName = ['/vols/Scratch/abaram/NN2020_BIDS/rawdata/' subjects{iSub} '/fmap/' subjects{iSub} '_phasediff.json'];
    savejson('',objPhase,opt);
end

%% dataset description

obj.Name             = 'Relational Structures in Reinforcement Learning';
obj.BIDSVersion      = '1.0.3';
obj.Authors      = {'Alon B Baram','Timothy H Muller','Hamed Nili','Mona Garvert','Timothy E J Behrens'};
obj.License          = 'CC0';
obj.Acknowledgements = 'We thank Jacob Bakermans, Anna Shpektor, Philipp Schwartenbeck, Avital Hahamy  and Shirley Mark for useful comments on earlier versions of the manuscript.';
opt.FileName = ['/vols/Scratch/abaram/NN2020_BIDS/rawdata/dataset_description.json'];
savejson('',obj,opt);


