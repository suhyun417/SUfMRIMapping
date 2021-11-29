% genFig_multipleFP_Revision_corrMap.m
%
% 2021/11/18 SHP
% Compute SU fMRI maps for each movie, then compare the results
%   - this is to reply to the reviewer's point during the revision for
%   Science Advances
% Written in Biowulf envioronment

clear all;

%% Settings
flagBiowulf = 1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        directory.projects = '/Volumes/PROJECTS';
        directory.procdata = '/Volumes/PROCDATA';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/Volumes/LIBRARY';
        addpath(fullfile(directory.library, 'matlab_utils'));
    else % on virtual machine
        directory.projects = '/nifvault/NIFVAULT/projects';
        directory.procdata = '/procdata';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/library';
        addpath(fullfile(directory.library, 'matlab_utils'));
    end
end

%% Load SU and fMRI tc
load(fullfile(directory.dataHome, 'matSDF_Movie123_allCells.mat'), 'matTS_FP') % matTS_FP.matNeuralRGR: [375×389 double]
load(fullfile(directory.dataHome, 'Art/fmritc4computeCorrcoeffCI_Art.mat')) % fmritc: time x movie x voxels

% get only face selective cells
load(fullfile(directory.dataHome, 'multipleFP_fsi.mat'))
locFaceCell =  find(fsi.matFSI(:,1)>0.33); % find(abs(fsi.matFSI(:,1))>0.33);

matNeuralRGR = matTS_FP.matNeuralRGR(:, locFaceCell);
% matTS_norm = zscore(matFR_TR); 
catAreaID = matTS_FP.catAreaID(locFaceCell);
catChanID = matTS_FP.catChanID(locFaceCell);
catSubjID = matTS_FP.catSubjID(locFaceCell);
setArea = matTS_FP.setArea;


%% Compute Correlation for each movie
tic;
corrMap_indMov = struct([]);

my_pool = parpool(4);
for iMovie = 1:3
    
    fprintf(1, '::::: Compute correlation for Movie %d :::::\n', iMovie); % Put some slider while corr is being calculated..
    
    curfmritc = squeeze(fmritc(:, iMovie, :)); % time x voxels
    
    matR_SU = NaN(size(curfmritc, 2), size(matNeuralRGR, 2));
    parfor iChan = 1:size(matNeuralRGR, 2)
        curneuralrgrs = matNeuralRGR(125*(iMovie-1)+1:125*iMovie, iChan); % time x neurons
        
        % Compute correlation
        fprintf(1, ':: chan #%d/%d  \n', iChan, size(matNeuralRGR, 2)); 
        
        [Rvals] = corr(curfmritc, curneuralrgrs, 'rows','complete', 'type', 'Spearman');
        
        matR_SU(:,iChan) = Rvals.*(-1); % because of MION
    end
    
    corrMap_indMov(iMovie).matR = matR_SU;
    
end
toc;
delete(my_pool)

% save(fullfile(directory.dataHome, 'multipleFP_Revision_corrMapIndMov.mat'), 'corrMap_indMov');


%% Compare fROI-versions
% Load the ROI indices
load(fullfile(directory.dataHome, 'Art/ROIs/Art_ROIs_set01_RH.mat')) %%s_ROIs_set00_RH.mat', nameSubjBOLD));

% % Load correlation matrix of all cortical cells
% load(sprintf('/procdata/parksh/_macaque/CorrMap_SU_AllCells%s_corticalFPMerged.mat', nameSubjBOLD), 'info*', 'corrMap_merged_FP'); %, 'info*', 'corrMap_Area', 'corrMap_merged');
% % load(sprintf('/procdata/parksh/_macaque/CorrMap_SU_AllCells%s_corticalFPMerged_pcares.mat', nameSubjBOLD), 'info*', 'corrMap_merged_FP'); %
% 
% dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

numROI = size(matROIIndices, 2);
nVoxPerROI = sum(matROIIndices);

for iMovie = 1:3
    matR_SU_all = corrMap_indMov(iMovie).matR; %corrMap_merged_FP.matR; %cat(2, corrMap_Area(setArea).matR); % except for no face patch neurons
    % catSubjID = corrMap_merged_FP.catSubjID; %cat(1, corrMap_Area(setArea).catSubjID);
    % catAreaID = corrMap_merged_FP.catAreaID; %floor(cat(1, corrMap_Area(setArea).catSubjID)./10);
    % catChanID = corrMap_merged_FP.catChanID;
    
    sumCorrROI = matR_SU_all'*matROIIndices;
    meanCorrROI = sumCorrROI./repmat(nVoxPerROI, size(sumCorrROI, 1), 1);
    
    corrMap_indMov(iMovie).meanCorrROI = meanCorrROI;
end


%% Compare SUfMRI maps from each movie and the entire set of three movies
load(fullfile(directory.dataHome, 'Art', 'Art_MaskArrays.mat'), 'movieDrivenAmp')
nx = 40; ny = 64; nz = 32;
nVox = nx*ny*nz;
brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
% matR_SU_all_brainmask = matR_SU_all(brainmask_vec,:); %matR_SU(brainmask_vec,:); % 27113 voxels

load(fullfile(directory.dataHome, 'Art', 'matR4clusteringmultiplepatches_faceselective.mat'))

setPair = nchoosek(1:3, 2);
for iPair = 1:length(setPair)
    [tempMatR] = corr(corrMap_indMov(setPair(iPair,1)).matR(brainmask_vec, :), ...
        corrMap_indMov(setPair(iPair,2)).matR(brainmask_vec, :),...
        'rows','complete', 'type', 'Spearman');
    setR = diag(tempMatR);
    ind = tril(true(size(tempMatR)), -1);
    setR_diffPairs = tempMatR(ind);
        
    corrMap_indMov_comparison(iPair).setMoviePair = setPair(iPair, :);
    corrMap_indMov_comparison(iPair).setR = setR;
    corrMap_indMov_comparison(iPair).setR_mean = mean(setR);    
    corrMap_indMov_comparison(iPair).setR_median = median(setR);  
    corrMap_indMov_comparison(iPair).setR_ste = std(setR)./sqrt(length(setR)-1);   
    corrMap_indMov_comparison(iPair).setR_diffPairs = setR_diffPairs;  
    
    % fROI-version
    [tempMatR_ROI] = corr(corrMap_indMov(setPair(iPair,1)).meanCorrROI', ...
        corrMap_indMov(setPair(iPair,2)).meanCorrROI',...
        'rows','complete', 'type', 'Spearman');
    ROI_setR = diag(tempMatR_ROI);
    ind = tril(true(size(tempMatR_ROI)), -1);
    ROI_setR_diffPairs = tempMatR_ROI(ind);

    corrMap_indMov_comparison(iPair).ROI_setR = ROI_setR;
    corrMap_indMov_comparison(iPair).ROI_setR_mean = mean(ROI_setR);    
    corrMap_indMov_comparison(iPair).ROI_setR_median = median(ROI_setR);  
    corrMap_indMov_comparison(iPair).ROI_setR_ste = std(ROI_setR)./sqrt(length(ROI_setR)-1);   
    corrMap_indMov_comparison(iPair).ROI_setR_diffPairs = ROI_setR_diffPairs; 
end    

% save(fullfile(directory.dataHome, 'multipleFP_Revision_corrMapIndMov.mat'), 'corrMap_indMov_comparison', '-append');

figure
for iP = 1:3
subplot(1,3,iP)
histogram(corrMap_indMov_comparison(iP).ROI_setR, 20)
end

%% Compare with the current 3-movie results
load(fullfile(directory.dataHome, 'Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat'), 'Clustering_meanROI')
load(fullfile(directory.dataHome, 'Art/Clustering_CorrMap_4FPs_faceselective_Movie123_probability.mat'), 'Clustering_brainmask', 'param*') 

% ROI sorting
setK = 2:20; %paramClustering_global.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS(:,iK) = sum(Clustering_brainmask.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
    matWSS_roi(:,iK) = sum(Clustering_meanROI.resultKMeans(iK).roi_sumD); % for grouping ROIs for visualization purposes
end

totalSS = Clustering_brainmask.totalSS_SU;
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

totalSS_roi = Clustering_meanROI.totalSS_roi;
propExplained_roi = (totalSS_roi-matWSS_roi)./totalSS_roi; 

% Get the ROI sorting
curK_roi = 9;
locMode_roi = find(propExplained_roi(:,curK_roi-1)==mode(propExplained_roi(:,curK_roi-1)));
locMin_roi = find(propExplained_roi(:,curK_roi-1)==min(propExplained_roi(:,curK_roi-1)));
[sortedClust_roi, indSortROI] = sort(Clustering_meanROI.resultKMeans(curK_roi-1).roi_indCluster(:, locMode_roi(1)));

% Plot
% Colormap
cMap_Area = [179 226 205; 141 160 203; 252 141 98; 231 41 138]./255; % 

fname = fullfile(directory.dataHome, 'Art/BCWYRColorMap.txt');
ttt = dlmread(fname);
cMap_corrSUMA = ttt(:, 1:3); % blue-cyan-white-yellow-red map for correlation
clear ttt

figure;
set(gcf, 'Color', 'w', 'Position', [550 200 860 790])
clf
for iM = 1:3
sp(iM) = subplot(1, 4, iM);
imagesc(corrMap_indMov(iM).meanCorrROI(:, indSortROI))
colormap(cMap_corrSUMA)
title(sprintf('Movie %d', iM))
end

sp(4) = subplot(1, 4, 4);
imagesc(Clustering_meanROI.matR(locFaceCell, indSortROI))
colormap(cMap_corrSUMA)
title('Three movies together')
set(sp, 'CLim', [-1 1].*0.7)

set(sp, 'YTick', [], 'XTick', []);

%% similarity of fROI-based correlation patterns
load(fullfile(directory.dataHome, 'multipleFP_Revision_corrMapIndMov.mat'), 'corrMap_indMov')
load(fullfile(directory.dataHome, 'Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat'), 'Clustering_meanROI')
% get only face selective cells
load(fullfile(directory.dataHome, 'multipleFP_fsi.mat'))
locFaceCell =  find(fsi.matFSI(:,1)>0.33); % find(abs(fsi.matFSI(:,1))>0.33);

for iM = 1:3
    
    % between each movie version and three-movie version
    [tempMatR_ROI] = corr(corrMap_indMov(iM).meanCorrROI', ...
        Clustering_meanROI.matR(locFaceCell, :)',...
        'rows','complete', 'type', 'Spearman');
    ROI_setR = diag(tempMatR_ROI);
    
    corrMap_indMov_comparison_withOriginal(iM).ROI_setR = ROI_setR;
    corrMap_indMov_comparison_withOriginal(iM).ROI_setR_mean = mean(ROI_setR);
    corrMap_indMov_comparison_withOriginal(iM).ROI_setR_median = median(ROI_setR);
    corrMap_indMov_comparison_withOriginal(iM).ROI_setR_ste = std(ROI_setR)./sqrt(length(ROI_setR)-1);
end

figure;
set(gcf, 'Position', [330 400 860 225])
for iM = 1:3
spp(iM) = subplot(1,3,iM);
histogram(corrMap_indMov_comparison_withOriginal(iM).ROI_setR, 20, 'EdgeColor', 'none'); hold on;
line([0 0], [0 30], 'LineStyle', '--', 'Color', 'k');
hold on;
plot(corrMap_indMov_comparison_withOriginal(iM).ROI_setR_median, 30, 'rv', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
% title(sprintf('Movie %d (5-min) and 15-min', iM))
end
set(spp,'XLim', [-1 1], 'Box', 'off', 'TickDir', 'out')
set(spp, 'YTick', 0:10:30, 'XColor', 'k', 'YColor', 'k', 'FontSize', 12)
% axis(spp, 'square')
print(gcf, fullfile(dirFig, 'multipleFP_Revision_figS_5minvs15min_withinCellbetweenMapCorr'), '-depsc')


figure;
set(gcf, 'Position', [330 400 860 530])
for iP = 1:3
spp(iP) = subplot(2,3,iP);
histogram(corrMap_indMov_comparison(iP).ROI_setR, 20, 'EdgeColor', 'none'); hold on;
line([0 0], [0 20], 'LineStyle', '--', 'Color', 'k');
hold on;
plot(corrMap_indMov_comparison(iP).ROI_setR_median, 20, 'rv', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
title(sprintf('Movie %d and %d', ...
    corrMap_indMov_comparison(iP).setMoviePair(1), corrMap_indMov_comparison(iP).setMoviePair(2)))
end
for iM = 1:3
spp(iM+3) = subplot(2,3,iM+3);
histogram(corrMap_indMov_comparison_withOriginal(iM).ROI_setR, 20, 'EdgeColor', 'none'); hold on;
line([0 0], [0 30], 'LineStyle', '--', 'Color', 'k');
hold on;
plot(corrMap_indMov_comparison_withOriginal(iM).ROI_setR_median, 20, 'rv', 'MarkerSize', 6, 'MarkerFaceColor', 'r')
title(sprintf('Movie %d and Original', iM))
end

set(spp,'XLim', [-1 1], 'Box', 'off', 'TickDir', 'out')


%% Heterogenity across cells?
mean(cat(2, corrMap_indMov_comparison.setR_diffPairs))

ans =

    0.0114    0.0322    0.0235

median(cat(2, corrMap_indMov_comparison.setR_diffPairs))

ans =

    0.0096    0.0369    0.0328

% across cell - 15min
tempR = corr(Clustering_meanROI.matR(locFaceCell, :)', 'rows', 'complete', 'type', 'Spearman');
ind = tril(true(size(tempR)), -1);
setR_diffPairs_org = tempR(ind);
median(setR_diffPairs_org)

ans =

    0.1344

mean(setR_diffPairs_org)

ans =

    0.1127

% across cell - 5 min
for iM = 1:3
tempR = corr(corrMap_indMov(iM).meanCorrROI', 'rows', 'complete', 'type', 'Spearman');
ind = tril(true(size(tempR)), -1);
setR_diffPairs = tempR(ind);

acrossCell(iM).setR_diffPairs = setR_diffPairs;
acrossCell(iM).medianR = median(setR_diffPairs);
acrossCell(iM).meanR = mean(setR_diffPairs);
end
    
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1700 390])
% imagesc(meanCorrROI')
% set(gca, 'XTick', cat(1, find(diff(catAreaID)>0), size(meanCorrROI, 1)), 'XTickLabel', corrMap_merged_FP.setArea)
% set(gca, 'YTick', 1:numROI, 'YTickLabel', paramROI.nameROI)
% set(gca, 'CLim', [-1 1].*0.5)
% set(gca, 'TickDir', 'out', 'Box', 'off')
% colorbar;
% title('mean correlation value for each ROI')
% 
% for iRoi = 1:numROI
%     setCorr = [];
%     setCorr = matR_SU_all(matROIIndices(:, iRoi)>0, :);
%     [y, i] = max(abs(setCorr));
%     ind = sub2ind(size(setCorr), i, 1:size(setCorr,2));
%     
%     corrROI(iRoi).setCorr = setCorr;
%     corrROI(iRoi).mean = mean(setCorr);
%     corrROI(iRoi).max = max(setCorr);
%     corrROI(iRoi).absmax = setCorr(ind);    
% end

