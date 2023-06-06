% exampleScript_SUCorrMap.m
%
% 2023/06/05 SHP
% - load the correlation matrix between individual neurons and fMRI ROIs
% - load the cluster identity of cells
% - example lines to reproduce Park et al. (2022) Fig. S8 (D) 

clear all;

% dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

load('/nifvault/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat', 'Clustering_meanROI')
load('/nifvault/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_probability.mat', 'Clustering_brainmask', 'param*')

%%
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

fROIs.nameROI = Clustering_meanROI.nameROI;
fROIs.orderROI_figS8 = indSortROI;

%% Fig 3B: correlation matrix before & after cell clustering

% Get the clustering results from voxel-based clustering
curK = 10; %13; % 9; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

% Colormap
cMap_Area = [179 226 205; 141 160 203; 252 141 98; 231 41 138]./255; % 
% cMap_Area = [91 148 203; 237 28 35; 248 148 29; 6 177 102]./255; % from Kenji's schematic
% cMap_Area(4, :) = cMap_Area(4, :).*0.7; % make the green a bit darker
% % orderArea = [4 1 2 3]; %ML-AF-AM-AAM
% % cMap_Area_MLfirst = cMap_Area(orderArea, :);

fname = '/nifvault/procdata/parksh/_macaque/Art/Anatomy/_suma/BCWYRColorMap.txt';
ttt = dlmread(fname);
cMap_corrSUMA = ttt(:, 1:3); % blue-cyan-white-yellow-red map for correlation
clear ttt

numROI = length(Clustering_meanROI.nameROI);
% orderROI = [1 2 22 3 4 35 34 12 13 14 29 30 6 7 8 36 9 10 11 32 15 5 23 26 37 27 28 16 17 33 18 19 20 21 31 24 25]; % 1:37;

tempS = struct([]);
% ordering of cells: ML-AF-AM-AAM, while maintaining the grouping
for iK = 1:curK
    tempS(iK).indSortChan_org = indSortChan(sortedClust==iK);
    tempS(iK).sortedClust_org = sortedClust(sortedClust==iK);
    tempS(iK).numCells = sum(sortedClust==iK);
    locML = find(Clustering_brainmask.infoCells.catAreaID(tempS(iK).indSortChan_org)==4);
    if isempty(locML)
        tempS(iK).indSortChan_reorder_MLfirst = tempS(iK).indSortChan_org;
    else
        tempS(iK).indSortChan_reorder_MLfirst = cat(1, tempS(iK).indSortChan_org(locML), tempS(iK).indSortChan_org(1:locML(1)-1));
    end
end
reorderCluster = [8 1 9 5 3 7 4 10 6 2]; % ORder cell groups based on the number of neurons in each group (from largest to smallest)
indSortChan_reorder = cat(1, tempS(reorderCluster).indSortChan_org);
indSortChan_reorder_MLfirst = cat(1, tempS(reorderCluster).indSortChan_reorder_MLfirst);
sortedClust_reorder = cat(1, tempS(reorderCluster).sortedClust_org);

%% panel D
fig3b2 = figure;
set(fig3b2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 350 1000]);

sp1 = subplot('Position', [0.08 0.05 0.62 0.9]);
imagesc(Clustering_meanROI.matR(indSortChan_reorder, indSortROI))
set(gca, 'CLim', [-1 1].*0.5)
locDiff = find(abs(diff(sortedClust_reorder))>0);
% locDiff_area = find(abs(diff(areaID_rand))>0);
set(sp1, 'YTick', [], 'YTickLabel', []) %locDiff+0.5, 'YTickLabel', [])
set(sp1, 'XTick', []) 
colormap(sp1, cMap_corrSUMA)
set(sp1,  'YColor', 'none', 'Box', 'off', 'XColor', 'none')
L = line(repmat([-0.8 37.5]', 1, length(locDiff)), [locDiff+0.5 locDiff+0.5]', 'Color', 'k',...
    'LineWidth', 2, 'LineStyle', ':');
set(sp1, 'XLim', [-0.8 37.5]);

sp2 = subplot('Position', [0.7 0.05 0.1 0.9]);
imagesc(Clustering_brainmask.infoCells.catAreaID(indSortChan_reorder)) %(indSortChan)')
set(sp2, 'XColor', 'none') %1, 'YTickLabel', 'Area info for each cell')
set(sp2, 'YTick', locDiff+0.5, 'YTickLabel', [])% line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k')
% xlabel('Number of cells in each cluster')
colormap(sp2, cMap_Area)
% spp2 = axes('Position', sp2.Position, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
% spp2.Clipping = 'off';
% spp2.XLim = sp2.XLim;
% spp2.YLim = sp2.YLim;
% spp2.YDir = sp2.YDir;
set(sp2,  'YColor', 'none', 'Box', 'off', 'XColor', 'none', 'XTick', [])
L2 = line(repmat([0.5;1.5], 1, length(locDiff)), [locDiff+0.5 locDiff+0.5]', 'Color', 'k',...
    'LineWidth', 2, 'LineStyle', ':');
set(sp2, 'XLim', [0.5 1.5]);


