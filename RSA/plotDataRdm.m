function plotDataRdm(rootData,subjects,glm,distType,analysisName,roi,fwhm)

% roi is a string with the vertex number and hemisphere, e.g. '80455lh'
% note the vertex number is expected in 0-based Freesufer indexes.

sl = 'surf'; % only data from cortical surfaces
modelRdm = getModelRdm(analysisName,false);

% initialise data RDM, excluding NaNs (this is so that for the relational
% structure analysis the RDM will be 8x8, not 12x12).

lengthRdmVec = length(modelRdm(~isnan(modelRdm))); % this is the length of the flattened RDM, after removing NaNs.
dataRdm = nan(length(subjects),lengthRdmVec); %
rdmInds = find(~isnan(modelRdm)); % non NaN elements of RDM. 
for iSub=1:length(subjects)
    sub=subjects{iSub};
    for iElem=1:length(rdmInds)
        dataFile = fullfile(rootData,'RDMs',glm,sub,'rdmElements',distType,sl,['elem' num2str(rdmInds(iElem)) '_smth' int2str(fwhm) '_' roi(end-1:end) '.mgh']);
        elementCurv = load_mgh(dataFile);
        dataRdm(iSub,iElem) = elementCurv(str2num(roi(1:end-2)) +1 ); % add 1 for conversion from 0-based freesurfer to 1-based Matlab
    end
end

dataRdmMean = mean(dataRdm,1);

%% plot data RDM

if strcmp(distType,'xRun')
    if strcmp(analysisName,'relationalStructure')
        lengthDiag = 8;
    elseif strcmp(analysisName,'visualIdentity')
        lengthDiag = 12;
    elseif strcmp(analysisName,'peXstructure')
        lengthDiag = 4;
    end
    
    figure;
    % off diagonal
    toPlot = squareform(dataRdmMean(1:lengthDiag*(lengthDiag-1)/2));
    
    % diagonal
    toPlot(logical(eye(lengthDiag))) = dataRdmMean(lengthDiag*(lengthDiag-1)/2+1:end);
    
elseif strcmp(distType,'corrWithinBlock')
    if strcmp('analysisName','latentVar_S-O')
        % first average the elements of S-O and O-S as this is
        % non-symmetric
        tmp = nan(3,1); 
        tmp(1) = (dataRdmMean(1) + dataRdmMean(3))/2;
        tmp(2) = (dataRdmMean(2) + dataRdmMean(5))/2;
        tmp(3) = (dataRdmMean(4) + dataRdmMean(6))/2;
        dataRdmMean = tmp;
    end
    lengthDiag = 3;

    figure;
    % only off diagonals
    toPlot = squareform(dataRdmMean); 
    % set diagonal to minimal value - this is only for visualisation
    % purposes as the diagonal is 0 by definition and i not used in the
    % stats.
    toPlot(logical(eye(lengthDiag))) = min(toPlot(toPlot>0));
    
end
if any(isnan(toPlot))
    % if there are NaNs, plot them in cian
    imagescwithnan(toPlot,hot,[0 1 1])
else
    % otherwise
    pcolor([toPlot,ones(size(toPlot,1),1).*mean(mean(toPlot)) ; ones(1,size(toPlot,2)+1).*mean(mean(toPlot))])
    colormap hot
    colorbar
end
axis image
axis ij
set(gca,'xtick',[])
set(gca,'xticklabel',[])
set(gca,'yticklabel',[])
set(gca,'ticklength',[0,0])
set(gca, 'FontSize', 18);