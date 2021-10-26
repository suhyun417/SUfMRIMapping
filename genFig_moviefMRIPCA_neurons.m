% genFig_moviefMRIPCA_neurons.m
%
% 2021/07/09 SHP: looking into face patch neurons correlated with fMRI
% principal components

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

%% Load data
load('/procdata/parksh/_macaque/ArtAva_movieTS_fMRI_concatTS_pca.mat', 'S'); % fMRI TS & PCs
load('/procdata/parksh/_macaque/matSDF_Movie123_allCells.mat', 'infoTS_subj', 'matTS_FP'); % neuronal time courses

% face selectivity
load('/procdata/parksh/_macaque/multipleFP_fsi.mat') 


%% re-form principal components
for iPC = 1:3
    for iS = 1:2
        
        tt = reshape(S(iS).resultsPCA_concat_brainmask.coeff(:,iPC), 118, 3);
        tt = cat(1, NaN(7, 3), tt);
        
        PC{iPC}(:,iS) = tt(:); %S(iS).resultsPCA_concat_brainmask.coeff(:,iPC);
    end
end

figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iPC = 1:3
SP(iPC) = subplot(3,1,iPC);
P{iPC} = plot(2.4:2.4:900, PC{iPC}, '-', 'LineWidth', 2);
axis tight
end

%% Correlation between PCs and neurons
matTS = matTS_FP.matNeuralRGR;
matTS_norm = zscore(matTS); 
catAreaID = matTS_FP.catAreaID;
catChanID = matTS_FP.catChanID;

[matR] = corr(PC{1}(:,1), matTS, 'type', 'Spearman', 'rows', 'complete'); 

figure;
set(gcf, 'Color', 'w')
histogram(matR, 50);
set(gca, 'XLim', [-1 1].*max(abs(get(gca, 'XLim'))))
xlabel('Spearman Correlation')
ylabel('Frequency')
title('Correlation between 1st PC of fMRI responses and face patch neurons')


[sortedR, indSortChan] = sort(matR, 'descend');

figure;
set(gcf, 'Color', 'w')
plot(sortedR, 'o-')
xlabel('Cells')
ylabel('Spearman Correlation')
title('Correlation between 1st PC of fMRI responses and face patch neurons: sorted')


















