% genFig_multipleFP_fig3fig4A_facecells.m
% 
% 2021/06/21 SHP
% Using only face cells
% 2020/09 SHP
% Multiple Face patch partially overlapping population manuscript
%  making figure 3: Clustering results using within-brain voxels. area
%  decomposition, etc.
%  plotting figure 4A: Cluster decomposition for each area

%% Clustering results using all the voxels, not ROI
clear all;

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat', 'Clustering_meanROI')
load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_faceselective_Movie123_probability.mat', 'Clustering_brainmask', 'param*') 

% get only face selective cells
load('/procdata/parksh/_macaque/multipleFP_fsi.mat')
locFaceCell =  find(fsi.matFSI(:,1)>0.33); % find(abs(fsi.matFSI(:,1))>0.33);

matR = Clustering_meanROI.matR(locFaceCell, :);
% matFR_TR = matTS_FP.matFR_TR(:, locFaceCell);
% matTS_norm = zscore(matFR_TR); 
% catAreaID = matTS_FP.catAreaID(locFaceCell);
% catChanID = matTS_FP.catChanID(locFaceCell);

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

% %% % Hierarchical clustering of ROIs % %%
% Y = pdist(Clustering_meanROI.matR');
% Z = linkage(Y, 'average');
% figure;
% [H, T] = dendrogram(Z, 0, 'Orientation', 'right');
% indY = str2num(get(gca, 'YTickLabel'));
% set(gca, 'YTickLabel', Clustering_meanROI.nameROI(indY))
% 
% indY = [7
%      9
%     35
%      4
%     13
%     12
%     34
%      1
%      2
%      3
%     17
%     33
%      5
%      8
%      6
%     16
%     21
%     23
%     36
%     37
%     26
%     22
%     10
%     11
%     18
%     19
%     20
%     29
%     32
%     14
%     15
%     24
%     25
%     27
%     28
%     31
%     30];
% 
% 
% 
% % figure;
% % set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% % plot(setK, propExplained'.*100, 'ko-', 'MarkerFaceColor', 'w'); hold on
% % xlabel('Number of cluster (K)')
% % ylabel('Explained variance (%)')
% % title('Clustering using all the voxels within brain')
% % set(gca, 'XTick', setK)
% % set(gca, 'TickDir', 'out', 'Box', 'off')
% % 
% % fig3b=figure;
% % set(fig3b, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 340 340])
% % curK=10;
% % plot(setK, propExplained.*100, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8); hold on
% % plot(curK, propExplained(:, curK-1).*100, 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 8)
% % xlim([2 20])
% % ylim([35 75])
% % set(gca, 'TickDir', 'out', 'LineWidth', 2, 'Box', 'off', 'TickLength', [.025 .05])
% % set(gca, 'YTick', 35:10:75)
% % 
% % % save
% % print(fig3b, fullfile(dirFig, 'expVar_KMeansClustering_brainmaskVoxels_square'), '-depsc')
% % 
% % 
% % figure;
% % set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% % plot(setK(1:end-1), diff(propExplained').*100)
% % hold on
% % plot(setK(1:end-1), mean(diff(propExplained').*100, 2), 'ko-', 'LineWidth', 2)
% % title('difference of explained variance for each K: using all the voxels within brain')

%% Fig 3B: correlation matrix before & after cell clustering

% Get the clustering results from voxel-based clustering
curK = 6; %8; %6; %10; %13; % 9; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

% Colormap
cMap_Area = [179 226 205; 141 160 203; 252 141 98; 231 41 138]./255; % 
% cMap_Area = [91 148 203; 237 28 35; 248 148 29; 6 177 102]./255; % from Kenji's schematic
% cMap_Area(4, :) = cMap_Area(4, :).*0.7; % make the green a bit darker
% % orderArea = [4 1 2 3]; %ML-AF-AM-AAM
% % cMap_Area_MLfirst = cMap_Area(orderArea, :);

fname = '/procdata/parksh/_macaque/Art/Anatomy/_suma/BCWYRColorMap.txt';
ttt = dlmread(fname);
cMap_corrSUMA = ttt(:, 1:3); % blue-cyan-white-yellow-red map for correlation
clear ttt

numROI = length(Clustering_meanROI.nameROI);
% orderROI = [1 2 22 3 4 35 34 12 13 14 29 30 6 7 8 36 9 10 11 32 15 5 23 26 37 27 28 16 17 33 18 19 20 21 31 24 25]; % 1:37;

% ordering of cells: ML-AF-AM-AAM, while maintaining the grouping
clear tempS
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
[~,reorderCluster] = sort(cat(1, tempS.numCells), 'descend');
% reorderCluster = 1:curK; %[8 1 9 5 3 7 4 10 6 2]; % ORder cell groups based on the number of neurons in each group (from largest to smallest)
indSortChan_reorder = cat(1, tempS(reorderCluster).indSortChan_org);
indSortChan_reorder_MLfirst = cat(1, tempS(reorderCluster).indSortChan_reorder_MLfirst);
sortedClust_reorder = cat(1, tempS(reorderCluster).sortedClust_org);

%% Fig 3B_1: Correlation matrix before clustering
for iArea = 1:4
    ttt(iArea).indChanArea = find((Clustering_brainmask.infoCells.catAreaID==iArea)>0);
    ttt(iArea).indChanArea_rand = ttt(iArea).indChanArea(randperm(length(ttt(iArea).indChanArea)));
end
indChanArea_rand = cat(1,  ttt([2 1 3 4]).indChanArea_rand);
areaID_rand = Clustering_brainmask.infoCells.catAreaID(indChanArea_rand);

fig3b_1 = figure;
set(fig3b_1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 350 1000]);

sp1 = subplot('Position', [0.08 0.05 0.62 0.9]);
imagesc(matR(indChanArea_rand, indSortROI))
set(gca, 'CLim', [-1 1].*0.5)
locDiff_area = find(abs(diff(areaID_rand))>0);
set(sp1, 'YTick', [], 'YTickLabel', []) %locDiff+0.5, 'YTickLabel', [])
set(sp1, 'XTick', []) 
colormap(sp1, cMap_corrSUMA)
set(sp1,  'YColor', 'none', 'Box', 'off', 'XColor', 'none')
L = line(repmat([-0.8 37.5]', 1, length(locDiff_area)), [locDiff_area+0.5 locDiff_area+0.5]', 'Color', 'k',...
    'LineWidth', 2, 'LineStyle', ':');
set(sp1, 'XLim', [-0.8 37.5]);

sp2 = subplot('Position', [0.7 0.05 0.1 0.9]);
imagesc(areaID_rand) %(indSortChan)')
set(sp2, 'XColor', 'none') %1, 'YTickLabel', 'Area info for each cell')
set(sp2, 'YTick', locDiff_area+0.5, 'YTickLabel', [])% line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k')
% xlabel('Number of cells in each cluster')
colormap(sp2, cMap_Area)
% spp2 = axes('Position', sp2.Position, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
% spp2.Clipping = 'off';
% spp2.XLim = sp2.XLim;
% spp2.YLim = sp2.YLim;
% spp2.YDir = sp2.YDir;
set(sp2,  'YColor', 'none', 'Box', 'off', 'XColor', 'none', 'XTick', [])
L2 = line(repmat([0.5;1.5], 1, length(locDiff_area)), [locDiff_area+0.5 locDiff_area+0.5]', 'Color', 'k',...
    'LineWidth', 2, 'LineStyle', ':');
set(sp2, 'XLim', [0.5 1.5]);

sp3 = subplot('Position', [0.8 0.05 0.06 0.9]);
% imagesc(setFSI(indSortChan_reorder))
% colormap(sp3, cMap_fsi)
% set(sp3, 'XColor', 'none') %1, 'YTickLabel', 'Area info for each cell')
set(sp3,  'YDir', 'reverse', 'YColor', 'none', 'Box', 'off', 'XColor', 'none', 'XTick', [])
L3 = line(repmat([0.5; 1], 1, length(locDiff_area)), [locDiff_area+0.5 locDiff_area+0.5]', 'Color', 'k',...
    'LineWidth', 2, 'LineStyle', ':');
set(sp3, 'YLim', [0.5 size(matR, 1)+0.5]);

% print(fig3b_1, fullfile(dirFig, 'Fig3B1_matR_meanROIbyCells_sortedArea_newColor_facecellsonly'), '-r200', '-dtiff');
 

%% Fig 3B_2: Correlation matrix after clustering
fig3b2 = figure;
set(fig3b2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 350 1000]);

sp1 = subplot('Position', [0.08 0.05 0.62 0.9]);
imagesc(matR(indSortChan_reorder, indSortROI))
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

sp3 = subplot('Position', [0.8 0.05 0.1 0.9]);
% imagesc(setFSI(indSortChan_reorder))
% colormap(sp3, cMap_fsi)
% set(sp3, 'XColor', 'none') %1, 'YTickLabel', 'Area info for each cell')
set(sp3,  'YDir', 'reverse', 'YColor', 'none', 'Box', 'off', 'XColor', 'none', 'XTick', [])
L3 = line(repmat([0.5; 1], 1, length(locDiff)), [locDiff+0.5 locDiff+0.5]', 'Color', 'k',...
    'LineWidth', 2, 'LineStyle', ':');
set(sp3, 'YLim', [0.5 size(matR, 1)+0.5]);

print(fig3b2, fullfile(dirFig, sprintf('Fig3B2_matR_meanROIbyCells_sortedGroup_K%d_newColor_facecells', curK)), '-r200', '-dtiff');

% % Mark fig 2 example neurons?
% setExampleCellIDs = {'33Dav', '130AFMoc', '097aMat', '10Dan'; ...
%     '27Dav', '065aTor', '39AMWas', '117AMMoc'; ...
%     '25Dav', '022bSpi', '51AMWas', '109AMMoc'; ...
%     '16Dav', '045aSpi', '33AMWas', '05Dan'; ...
%     '06Dav', '122AFMoc', '06AMWas', '115AMMoc'};
% 
% chanID_reorder = Clustering_meanROI.catChanID(indSortChan_reorder);
% locExampleCell = find(contains(Clustering_meanROI.catChanID(indSortChan_reorder), setExampleCellIDs')>0);
% 
% fig_temp = figure;
% set(fig_temp, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 350 1000]);
% sp1 = subplot('Position', [0.08 0.05 0.62 0.9]);
% imagesc(Clustering_meanROI.matR(indSortChan_reorder, indSortROI))
% set(sp1, 'YAxisLocation', 'right')
% set(gca, 'CLim', [-1 1].*0.5)
% locDiff = find(abs(diff(sortedClust_reorder))>0);
% L = line(repmat([-0.8 37.5]', 1, length(locDiff)), [locDiff+0.5 locDiff+0.5]', 'Color', 'k');
% colormap(sp1, cMap_corrSUMA)
% set(sp1, 'YTick', locExampleCell, 'YTIckLabel', chanID_reorder(locExampleCell))
% set(sp1, 'TickDir', 'out')
% set(sp1, 'XColor', 'none')

% % %Temp: fsi
% % load('/procdata/parksh/_macaque/multipleFP_fsi.mat') 
% % 
% % setFSI = fsi.matFSI(:,2); % taking only FSI
% % setFSI(abs(fsi.matFSI(:,2))<1) = 0; % no fingerprinting data available
% % setFSI(fsi.matFSI(:,2)<0) = 0.5; % not face selective
% % 
% % locDiff = find(abs(diff(sortedClust_reorder))>0);
% % cMap_fsi = [1 1 1; 0.7 0.7 0.7; 0 0 0]; % no data = white, not face selective = gray, face selective = black
% % 
% % figure(fig3b2);
% % sp3 = subplot('Position', [0.78 0.05 0.05 0.9]);
% % imagesc(setFSI(indSortChan_reorder))
% % colormap(sp3, cMap_fsi)
% % set(sp3, 'XColor', 'none') %1, 'YTickLabel', 'Area info for each cell')
% % set(sp3,  'YColor', 'none', 'Box', 'off', 'XColor', 'none', 'XTick', [])
% % L3 = line(repmat([0.5; 2.5], 1, length(locDiff)), [locDiff+0.5 locDiff+0.5]', 'Color', 'k',...
% %     'LineWidth', 2, 'LineStyle', ':');
% % set(sp3, 'XLim', [0.5 2.5]);
% % 
% % % print(fig3b2, fullfile(dirFig, 'Fig3B2_matR_meanROIbyCells_sortedGroup_K10_purplegreen_fsi'), '-r200', '-dtiff');

%% Fig 3C: average correlation map of each cell cluster: generate colormap with color-code of correlation in each fROI
% For each cell cluster, average correlation for each fROI across cells
% make it to "reordered" cluster ID
matAvgRforROI_Cluster = NaN(numROI, curK); % for easier ROI-color mapping in AFNI/SUMA, keep the original ROI order
for iK = 1:curK
    idK = reorderCluster(iK);
    matAvgRforROI_Cluster(:, iK) = mean(matR(tempS(idK).indSortChan_org, :))'; 
%     matAvgRforROI_Cluster_median(:, iK) = median(Clustering_meanROI.matR(tempS(idK).indSortChan_org, :))'; 
end

% Set the caxis limit
cmin = -0.3; %-0.4;
cmax = 0.3; %0.4;

% Convert correlation value to a scaled index given the current color axis
matColorIndex = fix(((matAvgRforROI_Cluster - cmin)./(cmax-cmin)*256)+1);
matColorIndex(matColorIndex > 256) = 256; % saturating it to maximum if it's larger than cmax

% Each cluster, save the color lookup table for 37 fROI using BCWYR colormap
fname = '/procdata/parksh/_macaque/Art/Anatomy/_suma/BCWYRColorMap.txt';
ttt = dlmread(fname);
cMap_corrSUMA = ttt(:, 1:3); % blue-cyan-white-yellow-red map for correlation
clear ttt

for iK = 1:curK
    clut = [];
    clut = cat(1, [1 1 1], cMap_corrSUMA(matColorIndex(:, iK), :));
    ind = (0:1:length(clut)-1)';
    clut = cat(2, clut, ind);
    
    % write a text file
    climstr = strrep(sprintf('%0.1f', abs(cmax)), '.', 'p');
    dlmwrite(sprintf('colorLUT_multipleFPSUMapping_onlyfacecells_fROI%d_cLim%s_ReorderedCellGroup%02d.txt', ...
        length(clut)-1, climstr, iK), clut);
end

[s, m, mid] = copyfile('./colorLUT*.txt', '/procdata/parksh/_macaque/Art/Anatomy/_suma');

clut_test = cMap_corrSUMA(matColorIndex(:, iK), :);
T = table([1:37]', string(Clustering_meanROI.nameROI), clut_test, ones([37 1]))
writetable(T, 'testClut_IndexedLabeled.txt', 'Delimiter', 'space', 'WriteVariableNames', false)

%% Fig 4A (or Fig 3D): composition of cell groups in each recording site
% Get the clustering results from voxel-based clustering
curK = 6; %8; %6; %10; %13; % 9; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

% % reorderCluster = 1:curK; %[8 1 9 5 3 7 4 10 6 2]; % ORder cell groups based on the number of neurons in each group (from largest to smallest)
% indSortChan_reorder = cat(1, tempS(reorderCluster).indSortChan_org);
% indSortChan_reorder_MLfirst = cat(1, tempS(reorderCluster).indSortChan_reorder_MLfirst);
% sortedClust_reorder = cat(1, tempS(reorderCluster).sortedClust_org);

cMap_cluster = pink(curK);

cellCountCluster_Area = NaN(curK, length(Clustering_brainmask.infoCells.setArea));
for iArea = 1:length(Clustering_brainmask.infoCells.setArea)
    compArea = [];
    compArea = sortedClust(ismember(indSortChan, find(Clustering_brainmask.infoCells.catAreaID == iArea)));
    cellCountCluster_Area(:, iArea) = histc(compArea, 1:curK);
end
cellCountCluster_Area_prop = cellCountCluster_Area./repmat(sum(cellCountCluster_Area), curK, 1);

% ordering of cells: ML-AF-AM-AAM, while maintaining the grouping
clear tempS
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
[~,reorderCluster] = sort(cat(1, tempS.numCells), 'descend');
% reorderCluster = [8 1 9 5 3 7 4 10 6 2]; % ORder cell groups based on the number of neurons in each group (from largest to smallest)
orderArea = [4 1 2 3]; %ML-AF-pAM-aAM

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 200 1000 290]);
bb = barh(cellCountCluster_Area_prop(reorderCluster, orderArea)', 'stacked'); %barh(cellCountCluster_Area_prop(:, orderArea)', 'stacked');
set(gca, 'TickDir', 'out', 'Box', 'off');
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', Clustering_brainmask.infoCells.setArea(orderArea))
set(gca, 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XColor', 'none')

lumFactor = linspace(0.3, 1.8, curK); %10);
for iB = 1:curK %10
bb(iB).FaceColor = 'flat';
bb(iB).CData = cMap_cluster(iB,:); 
% bb(iB).CData = cMap_Area_MLfirst*lumFactor(iB);
end

% add lines across bars
cumYData = cumsum(cat(1, bb.YData));
for iB = 1:curK %10
    xC = [bb(iB).XData(1)+bb(iB).BarWidth/2, bb(iB).XData(2)-bb(iB).BarWidth/2, bb(iB).XData(2)+bb(iB).BarWidth/2, ...
        bb(iB).XData(3)-bb(iB).BarWidth/2, bb(iB).XData(3)+bb(iB).BarWidth/2, bb(iB).XData(4)-bb(iB).BarWidth/2];
    yC = cumYData(iB, [1 2 2 3 3 4]);
    hold on
    line(reshape(yC, 2, 3), reshape(xC, 2, 3), 'Color', 'k', 'LineWidth', 0.5)
end

line([0.9 1], [5 5], 'Color', 'k', 'LineWidth', 3)
set(gca, 'YColor', 'none')
print(gcf, fullfile(dirFig, sprintf('cellGroupEachArea_facecells_barh_connected_nolabel_K%d_GroupReordered', curK)), '-r200', '-dtiff')
print(gcf, fullfile(dirFig, sprintf('cellGroupEachArea_facecells_barh_connected_nolabel_K%d_GroupReordered', curK)), '-depsc')

% color legend
fig_colorCluster = figure;
set(fig_colorCluster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 300 20 200]);
ax = subplot('Position', [0 0 1 1]);
image([1:curK]');
colormap(pink(curK));
axis off
line(repmat(get(gca, 'XLim')', 1, curK+1), repmat(0.5:1:curK+0.5, 2, 1), 'Color', 'k')
set(gca, 'XTick', [], 'YTick', [])
box on
print(fig_colorCluster, fullfile(dirFig, sprintf('colormap_pink%d_cluster', curK)), '-depsc')
print(fig_colorCluster, fullfile(dirFig, sprintf('colormap_pink%d_cluster', curK)), '-r200', '-dtiff')

%%% Fig 4A: cluster decomposition for each recording site (With reordered cell groups)
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 200 1000 290]);
% bb = barh(cellCountCluster_Area_prop(:, orderArea)', 'stacked');
% set(gca, 'TickDir', 'out', 'Box', 'off');
% set(gca, 'YDir', 'reverse')
% set(gca, 'YTickLabel', Clustering_brainmask.infoCells.setArea(orderArea))
% set(gca, 'LineWidth', 2, 'FontSize', 15)
% set(gca, 'XColor', 'none')
% cMap_cluster = pink(size(cellCountCluster_Area_prop, 1));
% lumFactor = linspace(0.3, 1.8, 10);
% for iB = 1:10
% bb(iB).FaceColor = 'flat';
% bb(iB).CData = cMap_cluster(iB,:); 
% % bb(iB).CData = cMap_Area_MLfirst*lumFactor(iB);
% end
% 
% % add lines across bars
% cumYData = cumsum(cat(1, bb.YData));
% for iB = 1:10
%     xC = [bb(iB).XData(1)+bb(iB).BarWidth/2, bb(iB).XData(2)-bb(iB).BarWidth/2, bb(iB).XData(2)+bb(iB).BarWidth/2, ...
%         bb(iB).XData(3)-bb(iB).BarWidth/2, bb(iB).XData(3)+bb(iB).BarWidth/2, bb(iB).XData(4)-bb(iB).BarWidth/2];
%     yC = cumYData(iB, [1 2 2 3 3 4]);
%     hold on
%     line(reshape(yC, 2, 3), reshape(xC, 2, 3), 'Color', 'k', 'LineWidth', 0.5)
% end
% 
% line([0.9 1], [5 5], 'Color', 'k', 'LineWidth', 3)
% set(gca, 'YColor', 'none')
% % print(gcf, fullfile(dirFig, 'cellGroupEachArea_barh_connected_nolabel'), '-r200', '-dtiff')
% % print(gcf, fullfile(dirFig, 'cellGroupEachArea_barh_connected_nolabel'), '-depsc')
% 
% % color legend
% fig_colorCluster = figure;
% set(fig_colorCluster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 300 20 200]);
% ax = subplot('Position', [0 0 1 1]);
% image([1:10]');
% colormap(pink(10));
% axis off
% line(repmat(get(gca, 'XLim')', 1, 11), repmat(0.5:1:10.5, 2, 1), 'Color', 'k')
% set(gca, 'XTick', [], 'YTick', [])
% box on
% % print(fig_colorCluster, fullfile(dirFig, 'colormap_pink10_cluster'), '-depsc')
% % print(fig_colorCluster, fullfile(dirFig, 'colormap_pink10_cluster'), '-r200', '-dtiff')




%%
% Fig3: previous format
% % Fig 3A
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 390 1160]);
% 
% sp1 = subplot('Position', [0.1 0.05 0.6 0.9]);
% imagesc(Clustering_meanROI.matR(indChanArea_rand, indY)) %indSortROI))
% set(gca, 'CLim', [-1 1].*0.5)
% locDiff_area = cat(1, 0, find(abs(diff(areaID_rand))>0), length(areaID_rand));
% set(sp1, 'YTick', [], 'YTickLabel', []) %locDiff+0.5, 'YTickLabel', [])
% set(sp1, 'XTick', []) 
% colormap(sp1, cMap_corrSUMA)
% set(gca, 'TickDir', 'out', 'Box', 'on')
% line(repmat(get(gca, 'XLim')', 1, length(locDiff_area)), [locDiff_area+0.5 locDiff_area+0.5]', 'Color', 'k', 'LineWidth', 0.5)
% 
% sp2 = subplot('Position', [0.75 0.05 0.1 0.9]);
% imagesc(areaID_rand) %(indSortChan)')
% set(sp2, 'XColor', 'none') %1, 'YTickLabel', 'Area info for each cell')
% set(sp2, 'YTick', locDiff_area+0.5, 'YTickLabel', [])% line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k')
% % xlabel('Number of cells in each cluster')
% colormap(sp2, cMap_Area)
% spp2 = axes('Position', sp2.Position, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
% spp2.Clipping = 'off';
% spp2.XLim = sp2.XLim;
% spp2.YLim = sp2.YLim;
% spp2.YDir = sp2.YDir;
% line(repmat([-0.7;1.5], 1, length(locDiff_area)), [locDiff_area+0.5 locDiff_area+0.5]', 'Color', 'k')
%  
% 
% % Fig 3B_2: Correlation matrix after clustering
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 390 1160]);
% 
% sp1 = subplot('Position', [0.1 0.05 0.6 0.9]);
% imagesc(Clustering_meanROI.matR(indSortChan_reorder, indSortROI))
% set(gca, 'CLim', [-1 1].*0.5)
% locDiff = cat(1, 0, find(abs(diff(sortedClust_reorder))>0), length(sortedClust_reorder));
% set(sp1, 'YTick', [], 'YTickLabel', []) %locDiff+0.5, 'YTickLabel', [])
% set(sp1, 'XTick', []) 
% colormap(sp1, cMap_corrSUMA)
% set(gca, 'TickDir', 'out', 'Box', 'on')
% line(repmat(get(gca, 'XLim')', 1, length(locDiff)), [locDiff+0.5 locDiff+0.5]', 'Color', 'k', 'LineWidth', 0.5)
% 
% sp2 = subplot('Position', [0.75 0.05 0.1 0.9]);
% imagesc(Clustering_brainmask.infoCells.catAreaID(indSortChan_reorder)) %(indSortChan)')
% set(sp2, 'XColor', 'none') %1, 'YTickLabel', 'Area info for each cell')
% set(sp2, 'YTick', locDiff+0.5, 'YTickLabel', [])% line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k')
% % xlabel('Number of cells in each cluster')
% colormap(sp2, cMap_Area)
% spp2 = axes('Position', sp2.Position, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
% spp2.Clipping = 'off';
% spp2.XLim = sp2.XLim;
% spp2.YLim = sp2.YLim;
% spp2.YDir = sp2.YDir;
% line(repmat([-0.7;1.5], 1, length(locDiff)), [locDiff+0.5 locDiff+0.5]', 'Color', 'k')
% %



%%%%%%%%%%%%%%%%%
%%%% Below: Older version %%%%%
%%%%%%%%%%%%%%%%%%
% sp1 = subplot('Position', [0.2 0.2 0.7 0.7]);
% imagesc(Clustering_meanROI.matR(indSortChan_reorder, orderROI)') %(indSortChan, orderROI)')
% set(sp1, 'CLim', [-1 1].*0.5)
% set(sp1, 'YTick', 1:numROI, 'YTickLabel', Clustering_meanROI.nameROI(orderROI))
% locDiff = cat(1, 0, find(diff(sortedClust)>0), length(sortedClust));
% set(sp1, 'XTick', locDiff+0.5, 'XTickLabel', locDiff)
% % title(sprintf('Clustered cells from 4 FPs using mean corr for each ROI: K=%d', curK))
% % xlabel('Cumulative number of cells')
% colormap(sp1, cMap_corrSUMA)
% set(gca, 'TickDir', 'out', 'Box', 'off')
% line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k', 'LineWidth', 0.5)
% 
% set(sp1, 'XTickLabel', [])
% spp1 = axes('Position', sp1.Position, 'Color', 'none', 'XAxisLocation', 'top');
% spp1.XLim = sp1.XLim;
% set(spp1, 'XTick', locDiff+0.5, 'XTickLabel', locDiff, 'TickDir', 'out', 'XTickLabel', [], 'Box', 'off', 'YColor', 'none')
% 
% sp2 = subplot('Position', [0.2 0.08 0.7 0.05]);
% imagesc(Clustering_brainmask.infoCells.catAreaID(indSortChan_reorder)') %(indSortChan)')
% set(sp2, 'YColor', 'none') %1, 'YTickLabel', 'Area info for each cell')
% set(sp2, 'XTick', locDiff+0.5, 'XTickLabel', [], 'TickDir', 'out'); %cat(1, locDiff(1), diff(locDiff)))
% % line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k')
% % xlabel('Number of cells in each cluster')
% colormap(sp2, cMap_Area)
% spp2 = axes('Position', sp2.Position, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
% spp2.Clipping = 'off';
% spp2.XLim = sp2.XLim;
% spp2.YLim = sp2.YLim;
% line([locDiff+0.5 locDiff+0.5]', repmat([0.5;3], 1, length(locDiff)), 'Color', 'k')

% 
% %% Fig 3A & B
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 400 910 675])
% 
% sp1 = subplot('Position', [0.2 0.2 0.7 0.7]);
% imagesc(Clustering_meanROI.matR(indSortChan_reorder, orderROI)') %(indSortChan, orderROI)')
% set(sp1, 'CLim', [-1 1].*0.5)
% set(sp1, 'YTick', 1:numROI, 'YTickLabel', Clustering_meanROI.nameROI(orderROI))
% locDiff = cat(1, 0, find(diff(sortedClust)>0), length(sortedClust));
% set(sp1, 'XTick', locDiff+0.5, 'XTickLabel', locDiff)
% % title(sprintf('Clustered cells from 4 FPs using mean corr for each ROI: K=%d', curK))
% % xlabel('Cumulative number of cells')
% colormap(sp1, cMap_corrSUMA)
% set(gca, 'TickDir', 'out', 'Box', 'off')
% line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k', 'LineWidth', 0.5)
% 
% set(sp1, 'XTickLabel', [])
% spp1 = axes('Position', sp1.Position, 'Color', 'none', 'XAxisLocation', 'top');
% spp1.XLim = sp1.XLim;
% set(spp1, 'XTick', locDiff+0.5, 'XTickLabel', locDiff, 'TickDir', 'out', 'XTickLabel', [], 'Box', 'off', 'YColor', 'none')
% 
% sp2 = subplot('Position', [0.2 0.08 0.7 0.05]);
% imagesc(Clustering_brainmask.infoCells.catAreaID(indSortChan_reorder)') %(indSortChan)')
% set(sp2, 'YColor', 'none') %1, 'YTickLabel', 'Area info for each cell')
% set(sp2, 'XTick', locDiff+0.5, 'XTickLabel', [], 'TickDir', 'out'); %cat(1, locDiff(1), diff(locDiff)))
% % line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k')
% % xlabel('Number of cells in each cluster')
% colormap(sp2, cMap_Area)
% spp2 = axes('Position', sp2.Position, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
% spp2.Clipping = 'off';
% spp2.XLim = sp2.XLim;
% spp2.YLim = sp2.YLim;
% line([locDiff+0.5 locDiff+0.5]', repmat([0.5;3], 1, length(locDiff)), 'Color', 'k')
% % colorbar;
% 
% print(gcf, fullfile(dirFig, 'matR_meanROIbyCells_KMeansClustering_brainmaskVoxels_K10_connected_nolabel', '-r200', '-dtiff'));
% 
%% each ROI color code from AFNI i64 colormap
fid = fopen('/procdata/parksh/_macaque/Art/ROIs/i64_colorscale.pal');
A = fscanf(fid, '%s');
fclose(fid);
matRGB = sscanf(A(8:end), '#%2x%2x%2x', [3 inf])';
matRGB = flipud(unique(matRGB, 'rows', 'stable'));
% 
orderROI = indSortROI'; %[1 2 22 3 4 35 34 12 13 14 29 30 6 7 8 36 9 10 11 32 15 5 23 26 37 27 28 16 17 33 18 19 20 21 31 24 25]; % 1:37;
% 
fig_colorROI = figure;
set(fig_colorROI, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 300 20 740]);
ax = subplot('Position', [0 0 1 1]);
image(orderROI');
colormap(matRGB./255);
axis off
% % print(fig_colorROI, fullfile(dirFig, 'colormap_ROIset01_reorder'), '-depsc')
print(fig_colorROI, fullfile(dirFig, 'colormap_ROIset01_reorder_ROIclustering'), '-r200', '-dtiff')
% 
% 
% %% Fig 3B bottom panel bar graph
% cellCountCluster_Area = NaN(curK, length(Clustering_brainmask.infoCells.setArea));
% for iArea = 1:length(Clustering_brainmask.infoCells.setArea)
%     compArea = [];
%     compArea = sortedClust(ismember(indSortChan, find(Clustering_brainmask.infoCells.catAreaID == iArea)));
%     cellCountCluster_Area(:, iArea) = histc(compArea, 1:curK);
% end
% cellCountCluster_Area_prop = cellCountCluster_Area./repmat(sum(cellCountCluster_Area), curK, 1);
% 
% 
% orderArea = [4 1 2 3]; %ML-AF-AM-AAM
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1318 468 910 120])
% for iC = 1:10
% ss(iC) = subplot('Position', [0.2+0.073*(iC-1) 0.05 0.05 0.8]);
% B(iC) = bar(cellCountCluster_Area_prop(iC, orderArea));
% end
% set(ss(:), 'TickDir', 'out', 'TickLength', [0.05 0.025], 'Box', 'off', 'XTick', [], 'XColor', 'k', 'YColor', 'k')
% set(B(:), 'FaceColor', 'flat', 'CData', cMap_Area_MLfirst, 'EdgeColor', 'none')
% for iS = 1:10
%     yT = ss(iS).YTick;
%     ss(iS).YTick = [yT(1) yT(end)];
% %     set(ss(iS), 'YTick', get(ss(iS), 'YLim'))
% end
% 
% 
% %% Fig 4A: cluster decomposition for each recording site
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 200 1000 290]);
% bb = barh(cellCountCluster_Area_prop(:, orderArea)', 'stacked');
% set(gca, 'TickDir', 'out', 'Box', 'off');
% set(gca, 'YDir', 'reverse')
% set(gca, 'YTickLabel', Clustering_brainmask.infoCells.setArea(orderArea))
% set(gca, 'LineWidth', 2, 'FontSize', 15)
% set(gca, 'XColor', 'none')
% cMap_cluster = pink(size(cellCountCluster_Area_prop, 1));
% lumFactor = linspace(0.3, 1.8, 10);
% for iB = 1:10
% bb(iB).FaceColor = 'flat';
% bb(iB).CData = cMap_cluster(iB,:); 
% % bb(iB).CData = cMap_Area_MLfirst*lumFactor(iB);
% end
% 
% % add lines across bars
% cumYData = cumsum(cat(1, bb.YData));
% for iB = 1:10
%     xC = [bb(iB).XData(1)+bb(iB).BarWidth/2, bb(iB).XData(2)-bb(iB).BarWidth/2, bb(iB).XData(2)+bb(iB).BarWidth/2, ...
%         bb(iB).XData(3)-bb(iB).BarWidth/2, bb(iB).XData(3)+bb(iB).BarWidth/2, bb(iB).XData(4)-bb(iB).BarWidth/2];
%     yC = cumYData(iB, [1 2 2 3 3 4]);
%     hold on
%     line(reshape(yC, 2, 3), reshape(xC, 2, 3), 'Color', 'k', 'LineWidth', 0.5)
% end
% 
% line([0.9 1], [5 5], 'Color', 'k', 'LineWidth', 3)
% set(gca, 'YColor', 'none')
% % print(gcf, fullfile(dirFig, 'cellGroupEachArea_barh_connected_nolabel'), '-r200', '-dtiff')
% % print(gcf, fullfile(dirFig, 'cellGroupEachArea_barh_connected_nolabel'), '-depsc')
% 
% % color legend
% fig_colorCluster = figure;
% set(fig_colorCluster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 300 20 200]);
% ax = subplot('Position', [0 0 1 1]);
% image([1:10]');
% colormap(pink(10));
% axis off
% line(repmat(get(gca, 'XLim')', 1, 11), repmat(0.5:1:10.5, 2, 1), 'Color', 'k')
% set(gca, 'XTick', [], 'YTick', [])
% box on
% % print(fig_colorCluster, fullfile(dirFig, 'colormap_pink10_cluster'), '-depsc')
% % print(fig_colorCluster, fullfile(dirFig, 'colormap_pink10_cluster'), '-r200', '-dtiff')


