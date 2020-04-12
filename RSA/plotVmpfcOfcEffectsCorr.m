function plotVmpfcOfcEffectsCorr(rootData,subjects,fwhm)

effects = nan(length(subjects),2);
for iSub=1:length(subjects)
    sub=subjects{iSub};
    tmp = load_mgh(fullfile(rootData,'RDMs','GLM3',sub,'rsaContrasts','xRun','surf',...
    ['peXstructure_smth' int2str(fwhm) '_lh.mgh']));
    % look at the peak of vmpfc - vertex 110384
    effects(iSub,1) = tmp(110384 + 1); % +1 vor transition betwen Freesurfer and Matlab indeces. 
    tmp = load_mgh(fullfile(rootData,'RDMs','GLM2',sub,'rsaContrasts','corrWithinBlock','surf',...
    ['latentVar_all_smth' int2str(fwhm) '_lh.mgh']));
    % look at the peak of OFC - vertex 110371
    effects(iSub,2) = tmp(110371 + 1); % +1 vor transition betwen Freesurfer and Matlab indeces. 
end

% calc p value using permutation test
[rho_true, p_param] = corr(effects(:,1),effects(:,2));
nPerm = 10000;
rhos = nan(nPerm,1);
for iPerm=1:nPerm
    rhos(iPerm) = corr(effects(:,1),effects(randperm(28),2));
end
p_perm = sum(rho_true<rhos)/length(rho);

% plot
figure;
hold on
scatter(effects(:,1),effects(:,2),300,'.');
lsline
str=sprintf('r= %1.2f, p= %1.3f',rho_true(1,2),p_perm(1,2));
T = text(min(get(gca, 'xlim')), max(get(gca, 'ylim')), str);
set(T, 'fontsize', 18, 'verticalalignment', 'top', 'horizontalalignment', 'left');
set(gca,'fontsize', 18)
xlabel('vmPFC pe x structure')
ylabel('OFC latent variable')