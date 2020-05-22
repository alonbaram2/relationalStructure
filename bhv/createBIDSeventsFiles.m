clear
rootScripts = '/home/fs0/abaram/scripts/NN2020';
rootData    = '/home/fs0/abaram/scratch/NN2020';
subjects = {'sub-01','sub-02','sub-03','sub-04','sub-05','sub-06','sub-07','sub-11','sub-12','sub-13','sub-14','sub-15','sub-16','sub-17','sub-18','sub-19',...
    'sub-20','sub-21','sub-22','sub-23','sub-24','sub-25','sub-26','sub-27','sub-28'};

nSub = length(subjects);
nRun  = 2; % two independent runs, each with the 4 conditions
nCond = 4; % the 4 experimental conditions
nTrialsPerBlock = 30;

%% 

for iSub = 1: length(subjects)
    sub = subjects{iSub};
    f = [rootData '/bhv/raw/scanDay/' sub '/scan.mat'];
    load(f);
    for r = 1:nRun % two runs, each with the 4 conditions
        for c = 1:nCond % the 4 experimental conditions, already reordered:
                       % 1: stimSet 1, -Corr.
                       % 2: stimSet 2, -Corr. 
                       % 3: stimSet 1, +Corr. 
                       % 4: stimSet 2, +Corr.
                       
                       
                       % each trial has two events: stim and outcome
                       % presentations. So we need to duplicate the other
                       % fields
                       cue     = repelem(data.run{r}.cond{c}.cues,2);
                       outcome = repelem(data.run{r}.cond{c}.rewards,2);
                       choice  = repelem(data.run{r}.cond{c}.choices,2);
                       prob    = repelem(data.run{r}.cond{c}.probs,2);
                       reaction_time  = repelem(data.run{r}.cond{c}.RT,2);
                       trial  = repelem([1:nTrialsPerBlock]',2);
                                              
                       onsetStim       = data.run{r}.cond{c}.stimOn - data.run{r}.cond{c}.tTrig1;
                       onsetOutcome    = data.run{r}.cond{c}.feedbackOn - data.run{r}.cond{c}.tTrig1;
                       
                       onset = sort(vertcat(onsetStim,onsetOutcome));
                       
                       % durations are constant for all reials - take from
                       % the first trial
                       durationStim    = data.run{r}.cond{c}.choicePossibleOn(1) - data.run{r}.cond{c}.stimOn(1);                       
                       durationOutcome = data.run{r}.cond{c}.ITIfixationOn(1) - data.run{r}.cond{c}.feedbackOn(1);      
                       duration = repmat(vertcat(durationStim,durationOutcome),[nTrialsPerBlock,1]);
                       
                       % create indicator variable - is the event a stim or
                       % outcome presentation.
                       stimOrOutcome = repmat({'stim';'outcome'},[nTrialsPerBlock,1]);
                       
                       T = table(onset,duration,trial,stimOrOutcome,cue,outcome,prob,choice,reaction_time);
                       writetable(T,['/vols/Scratch/abaram/NN2020_BIDS/' subjects{iSub} ...
                           '/func/' subjects{iSub} '_task-cond0' num2str(c) '_run-0' num2str(r) '_events.tsv'],'FileType','text','Delimiter','\t')
        end
    end
end
