function vecRDM = getModelRdm(analysisName,plotFlag)

switch analysisName
    
    case 'visualIdentity'
        labels = {'A-','B-','C0','D-','E-','F0',...% Neg corr
            'A+','B+','C0','D+','E+','F0'}; % pos corr
        
        RDM = ones(12) - eye(12); % diag is 0 - same stims => 0 dist. this should give positive loadings for visual areas.
        for i=1:6
            RDM(i,i+6) = 0;  RDM(i+6,i) = 0;% same stims in different runs
        end
        diag = zeros(1,12);
        
    case 'relationalStructure'
        % same vs different structure, only stimuli A & B (separate regressor
        % for each). Data RDM also includes [C] trials, which need to be
        % ignored here.
        
        % namnes of different conditions. ABC and DEF are the two stimuli
        % sets.
        labels = {'A-','B-','C0','D-','E-','F0',...% Neg corr
            'A+','B+','C0','D+','E+','F0'}; % pos corr
        
        RDM = ones(12);
        %         % diag is zeros in order to make this a distance matrix that can be given to
        %         % squareform as input. (this will be ignored later, and we will set the
        %         % diag again)
        %         for i=1:length(RDM)
        %             RDM(i,i)=0;
        %         end
        % -corr to -corr
        RDM(1,2) = 0;
        RDM(4,5) = 0;
        RDM(1,4) = 0;
        RDM(1,5) = 0;
        RDM(2,4) = 0;
        RDM(2,5) = 0;
        % +corr to +corr
        RDM(7,8) = 0;
        RDM(10,11) = 0;
        RDM(7,10) = 0;
        RDM(7,11) = 0;
        RDM(8,10) = 0;
        RDM(8,11) = 0;
        
        % get rid of all unrelated (stim C)
        RDM(3,:) = nan;
        RDM(6,:) = nan;
        RDM(9,:) = nan;
        RDM(12,:) = nan;
        RDM(:,3) = nan;
        RDM(:,6) = nan;
        RDM(:,9) = nan;
        RDM(:,12) = nan;
        
        RDM = makeRdmSymmetric(RDM);
        
        % diagonal elements
        diag = [0,0,nan,0,0,nan,0,0,nan,0,0,nan];
        
    case 'relationalStructure_control1'
        % like the previous case but ignoring all same-stimulus elements.
        
        labels = {'A-','B-','C0','D-','E-','F0',...% Neg corr
            'A+','B+','C0','D+','E+','F0'}; % pos corr
        
        RDM = ones(12);
        %         % diag is zeros in order to make this a distance matrix that can be given to
        %         % squareform as input. (this will be ignored later, and we will set the
        %         % diag again)
        %         for i=1:length(RDM)
        %             RDM(i,i)=0;
        %         end
        
        % same stim, diff structure elements are ignored
        for i=1:length(RDM)/2
            RDM(i,i+length(RDM)/2)=nan;
        end
        % -corr to -corr
        RDM(1,2) = 0;
        RDM(4,5) = 0;
        RDM(1,4) = 0;
        RDM(1,5) = 0;
        RDM(2,4) = 0;
        RDM(2,5) = 0;
        % +corr to +corr
        RDM(7,8) = 0;
        RDM(10,11) = 0;
        RDM(7,10) = 0;
        RDM(7,11) = 0;
        RDM(8,10) = 0;
        RDM(8,11) = 0;
        
        % get rid of all unrelated (stim C)
        RDM(3,:) = nan;
        RDM(6,:) = nan;
        RDM(9,:) = nan;
        RDM(12,:) = nan;
        RDM(:,3) = nan;
        RDM(:,6) = nan;
        RDM(:,9) = nan;
        RDM(:,12) = nan;
        
        RDM = makeRdmSymmetric(RDM);
        
        % Ignore diagonal elements
        diag = nan(1,12);
    case 'relationalStructure_control2'
        % like the first case but ignoring all within stimuli-set elements.
        
        labels = {'A-','B-','C0','D-','E-','F0',...% Neg corr
            'A+','B+','C0','D+','E+','F0'}; % pos corr
        
        RDM = ones(12);
        
        % get rid of within stimuli-set elements
        
        RDM(1:3,1:3) = nan;
        RDM(4:6,4:6) = nan;
        RDM(7:9,7:9) = nan;
        RDM(10:12,10:12) = nan;
        RDM(1:3,7:9) = nan;
        RDM(4:6,10:12) = nan;
        
        
        %         % diag is zeros in order to make this a distance matrix that can be given to
        %         % squareform as input. (this will be ignored later, and we will set the
        %         % diag again)
        %         for i=1:length(RDM)
        %             RDM(i,i)=0;
        %         end
        
        % same stim, diff structure elements are ignored
        for i=1:length(RDM)/2
            RDM(i,i+length(RDM)/2)=nan;
        end
        % -corr to -corr across diff blocks (diff stim sets)
        RDM(1,4) = 0;
        RDM(1,5) = 0;
        RDM(2,4) = 0;
        RDM(2,5) = 0;
        % +corr to +corr across diff blocks (diff stim sets)
        RDM(7,10) = 0;
        RDM(7,11) = 0;
        RDM(8,10) = 0;
        RDM(8,11) = 0;
        
        % get rid of all unrelated (stim C)
        RDM(3,:) = nan;
        RDM(6,:) = nan;
        RDM(9,:) = nan;
        RDM(12,:) = nan;
        RDM(:,3) = nan;
        RDM(:,6) = nan;
        RDM(:,9) = nan;
        RDM(:,12) = nan;
        
        RDM = makeRdmSymmetric(RDM);
        
        % Ignore diagonal elements
        diag = nan(1,12);
        
    case 'peXstructure'
        labels = {'SS1_St1','SS2_St1','SS1_St2','SS2_St2'}; % SS1_St1 is stimuli set 1, structure 1
        RDM = nan(4);
        RDM(1,2) = 0; % -corr 2 -corr
        RDM(1,3) = 1; % -corr 2 +corr
        RDM(1,4) = 1; % -corr 2 +corr
        RDM(2,3) = 1; % -corr 2 +corr
        RDM(2,4) = 1; % -corr 2 +corr
        RDM(3,4) = 0; % +corr 2 +corr
        
        RDM = makeRdmSymmetric(RDM);
        diag = zeros(1,4); % all diagonal elements are same structure.
        
    case 'latentVar_all'
        labels = {'A_{S}','B_{S}','C_{S}','A_{O}','B_{O}','C_{O}'}; %
        RDM = nan(6);
        RDM(1,2) = 0; % (A_S,B_S)
        RDM(1,3) = 1; % (A_S,C_S)
        RDM(2,3) = 1; % (B_S,C_S)
        RDM(4,5) = 0; % (A_O,B_O)
        RDM(4,6) = 1; % (A_O,C_O)
        RDM(5,6) = 1; % (B_O,C_O)
        RDM(1,5) = 0; % (A_S,B_O)
        RDM(2,4) = 0; % (B_S,A_O)
        RDM(1,6) = 1; % (A_S,C_O)
        RDM(2,6) = 1; % (B_S,C_O)
        RDM(3,4) = 1; % (C_S,A_O)
        RDM(3,5) = 1; % (C_S,B_O)        
        RDM = makeRdmSymmetric(RDM);
        diag = []; % symmetric matrix, as we've used regular correlation ditance for the within-block analyses.
    case 'latentVar_S-S'
        labels = {'A_{S}','B_{S}','C_{S}','A_{O}','B_{O}','C_{O}'}; %
        RDM = nan(6);
        RDM(1,2) = 0; % (A_S,B_S)
        RDM(1,3) = 1; % (A_S,C_S)
        RDM(2,3) = 1; % (B_S,C_S)
        RDM = makeRdmSymmetric(RDM);        
        diag = []; % symmetric matrix, as we've used regular correlation ditance for the within-block analyses.
        
    case 'latentVar_O-O'
        labels = {'A_{S}','B_{S}','C_{S}','A_{O}','B_{O}','C_{O}'}; %
        RDM = nan(6);
        RDM(4,5) = 0; % (A_O,B_O)
        RDM(4,6) = 1; % (A_O,C_O)
        RDM(5,6) = 1; % (B_O,C_O)
        RDM = makeRdmSymmetric(RDM);        
        diag = []; % symmetric matrix, as we've used regular correlation ditance for the within-block analyses.
        
    case 'latentVar_S-O'
        labels = {'A_{S}','B_{S}','C_{S}','A_{O}','B_{O}','C_{O}'}; %
        RDM = nan(6);
        RDM(1,5) = 0; % (A_S,B_O)
        RDM(2,4) = 0; % (B_S,A_O)
        RDM(1,6) = 1; % (A_S,C_O)
        RDM(2,6) = 1; % (B_S,C_O)
        RDM(3,4) = 1; % (C_S,A_O)
        RDM(3,5) = 1; % (C_S,B_O)
        RDM = makeRdmSymmetric(RDM);        
        diag = []; % symmetric matrix, as we've used regular correlation ditance for the within-block analyses.
        
end

vecRDM = vertcat(vectoriseRDM(RDM),diag');


if plotFlag
    toPlot = RDM;
    if ~isempty(diag)
        toPlot(logical(eye(length(toPlot)))) = diag;
    end
    figure;
    imagescwithnan(toPlot,hot,[0 1 1])
    
    ax = gca();
    set(gca, 'FontSize', 18);
    ax.XTick=1:length(labels);
    ax.YTick=1:length(labels);
    ax.XTickLabels = labels;
    ax.YTickLabels = labels;
    ax.XTickLabelRotation = 45;
    
    axis square
    title(analysisName)
end

function RDM = makeRdmSymmetric(RDM)
% make symmetric
for i=1:length(RDM)
    for j=i:length(RDM)
        if i<j
            RDM(j,i)=RDM(i,j);
        end
    end
end

function v = vectoriseRDM(RDM)
if ~any(any(isnan(RDM))) % if there are nans, matlab will treat a symmetrical matrix as asymmetric!
    if any(any(RDM ~= RDM')) || size(RDM,1) ~= size(RDM,2)
        error('RDM must be square and symmetric')
    end
end
if any(diag(RDM))
    warning('non-zero elements on diagonal of RDM!')
end

v=RDM(tril(true(length(RDM)),-1));

