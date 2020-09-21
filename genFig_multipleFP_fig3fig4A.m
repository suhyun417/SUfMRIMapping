% genFig_multipleFP_fig3fig4A.m
% 
% 2020/09 SHP
% Multiple Face patch partially overlapping population manuscript
%  making figure 3: Clustering results using within-brain voxels. area
%  decomposition, etc.
%  plotting figure 4A: Cluster decomposition for each area

%% Clustering results using all the voxels, not ROI
clear all;

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat', 'Clustering_meanROI')
load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_probability.mat', 'Clustering_brainmask', 'param*') 


%%
setK = 2:20; %paramClustering_global.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS(:,iK) = sum(Clustering_brainmask.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
end

totalSS = Clustering_brainmask.totalSS_SU;
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% plot(setK, propExplained'.*100, 'ko-', 'MarkerFaceColor', 'w'); hold on
% xlabel('Number of cluster (K)')
% ylabel('Explained variance (%)')
% title('Clustering using all the voxels within brain')
% set(gca, 'XTick', setK)
% set(gca, 'TickDir', 'out', 'Box', 'off')
% 
% fig3b=figure;
% set(fig3b, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 340 340])
% curK=10;
% plot(setK, propExplained.*100, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8); hold on
% plot(curK, propExplained(:, curK-1).*100, 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 8)
% xlim([2 20])
% ylim([35 75])
% set(gca, 'TickDir', 'out', 'LineWidth', 2, 'Box', 'off', 'TickLength', [.025 .05])
% set(gca, 'YTick', 35:10:75)
% 
% % save
% print(fig3b, fullfile(dirFig, 'expVar_KMeansClustering_brainmaskVoxels_square'), '-depsc')
% 
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% plot(setK(1:end-1), diff(propExplained').*100)
% hold on
% plot(setK(1:end-1), mean(diff(propExplained').*100, 2), 'ko-', 'LineWidth', 2)
% title('difference of explained variance for each K: using all the voxels within brain')


%% Fig 3A & B: Correlation matrix
% Get the clustering results from voxel-based clustering
curK = 10; %13; % 9; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

% Colormap
cMap_Area = [91 148 203; 237 28 35; 248 148 29; 6 177 102]./255; % from Kenji's schematic
cMap_Area(4, :) = cMap_Area(4, :).*0.7; % make the green a bit darker
orderArea = [4 1 2 3]; %ML-AF-AM-AAM
cMap_Area_MLfirst = cMap_Area(orderArea, :);

fname = '/procdata/parksh/_macaque/Art/Anatomy/_suma/BCWYRColorMap.txt';
ttt = dlmread(fname);
cMap_corrSUMA = ttt(:, 1:3); % blue-cyan-white-yellow-red map for correlation
clear ttt

numROI = length(Clustering_meanROI.nameROI);
orderROI = [1 2 22 3 4 35 34 12 13 14 29 30 6 7 8 36 9 10 11 32 15 5 23 26 37 27 28 16 17 33 18 19 20 21 31 24 25]; % 1:37;

% ordering of cells: ML-AF-AM-AAM, while maintaining the grouping
for iK = 1:curK
    tempS(iK).indSortChan_org = indSortChan(sortedClust==iK);
    locML = find(Clustering_brainmask.infoCells.catAreaID(tempS(iK).indSortChan_org)==4);
    if isempty(locML)
        tempS(iK).indSortChan_reorder = tempS(iK).indSortChan_org;
    else
        tempS(iK).indSortChan_reorder = cat(1, tempS(iK).indSortChan_org(locML), tempS(iK).indSortChan_org(1:locML(1)-1));
    end
end
indSortChan_reorder = cat(1, tempS.indSortChan_reorder);

%% Fig 3A & B
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 400 910 675])

sp1 = subplot('Position', [0.2 0.2 0.7 0.7]);
imagesc(Clustering_meanROI.matR(indSortChan_reorder, orderROI)') %(indSortChan, orderROI)')
set(sp1, 'CLim', [-1 1].*0.5)
set(sp1, 'YTick', 1:numROI, 'YTickLabel', Clustering_meanROI.nameROI(orderROI))
locDiff = cat(1, 0, find(diff(sortedClust)>0), length(sortedClust));
set(sp1, 'XTick', locDiff+0.5, 'XTickLabel', locDiff)
% title(sprintf('Clustered cells from 4 FPs using mean corr for each ROI: K=%d', curK))
% xlabel('Cumulative number of cells')
colormap(sp1, cMap_corrSUMA)
set(gca, 'TickDir', 'out', 'Box', 'off')
line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k', 'LineWidth', 0.5)

set(sp1, 'XTickLabel', [])
spp1 = axes('Position', sp1.Position, 'Color', 'none', 'XAxisLocation', 'top');
spp1.XLim = sp1.XLim;
set(spp1, 'XTick', locDiff+0.5, 'XTickLabel', locDiff, 'TickDir', 'out', 'XTickLabel', [], 'Box', 'off', 'YColor', 'none')

sp2 = subplot('Position', [0.2 0.08 0.7 0.05]);
imagesc(Clustering_brainmask.infoCells.catAreaID(indSortChan_reorder)') %(indSortChan)')
set(sp2, 'YColor', 'none') %1, 'YTickLabel', 'Area info for each cell')
set(sp2, 'XTick', locDiff+0.5, 'XTickLabel', [], 'TickDir', 'out'); %cat(1, locDiff(1), diff(locDiff)))
% line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k')
% xlabel('Number of cells in each cluster')
colormap(sp2, cMap_Area)
spp2 = axes('Position', sp2.Position, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
spp2.Clipping = 'off';
spp2.XLim = sp2.XLim;
spp2.YLim = sp2.YLim;
line([locDiff+0.5 locDiff+0.5]', repmat([0.5;3], 1, length(locDiff)), 'Color', 'k')
% colorbar;

print(gcf, fullfile(dirFig, 'matR_meanROIbyCells_KMeansClustering_brainmaskVoxels_K10_connected_nolabel', '-r200', '-dtiff'));

%% each ROI color code from AFNI i64 colormap
fid = fopen('/procdata/parksh/_macaque/Art/ROIs/i64_colorscale.pal');
A = fscanf(fid, '%s');
fclose(fid);
matRGB = sscanf(A(8:end), '#%2x%2x%2x', [3 inf])';
matRGB = flipud(unique(matRGB, 'rows', 'stable'));

orderROI = [1 2 22 3 4 35 34 12 13 14 29 30 6 7 8 36 9 10 11 32 15 5 23 26 37 27 28 16 17 33 18 19 20 21 31 24 25]; % 1:37;

fig_colorROI = figure;
set(fig_colorROI, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 300 20 740]);
ax = subplot('Position', [0 0 1 1]);
image(orderROI');
colormap(matRGB./255);
axis off
% print(fig_colorROI, fullfile(dirFig, 'colormap_ROIset01_reorder'), '-depsc')
% print(fig_colorROI, fullfile(dirFig, 'colormap_ROIset01_reorder'), '-r200', '-dtiff')


%% Fig 3B bottom panel bar graph
cellCountCluster_Area = NaN(curK, length(Clustering_brainmask.infoCells.setArea));
for iArea = 1:length(Clustering_brainmask.infoCells.setArea)
    compArea = [];
    compArea = sortedClust(ismember(indSortChan, find(Clustering_brainmask.infoCells.catAreaID == iArea)));
    cellCountCluster_Area(:, iArea) = histc(compArea, 1:curK);
end
cellCountCluster_Area_prop = cellCountCluster_Area./repmat(sum(cellCountCluster_Area), curK, 1);


orderArea = [4 1 2 3]; %ML-AF-AM-AAM
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1318 468 910 120])
for iC = 1:10
ss(iC) = subplot('Position', [0.2+0.073*(iC-1) 0.05 0.05 0.8]);
B(iC) = bar(cellCountCluster_Area_prop(iC, orderArea));
end
set(ss(:), 'TickDir', 'out', 'TickLength', [0.05 0.025], 'Box', 'off', 'XTick', [], 'XColor', 'k', 'YColor', 'k')
set(B(:), 'FaceColor', 'flat', 'CData', cMap_Area_MLfirst, 'EdgeColor', 'none')
for iS = 1:10
    yT = ss(iS).YTick;
    ss(iS).YTick = [yT(1) yT(end)];
%     set(ss(iS), 'YTick', get(ss(iS), 'YLim'))
end


%% Fig 4A: cluster decomposition for each recording site
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 200 1000 290]);
bb = barh(cellCountCluster_Area_prop(:, orderArea)', 'stacked');
set(gca, 'TickDir', 'out', 'Box', 'off');
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', Clustering_brainmask.infoCells.setArea(orderArea))
set(gca, 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XColor', 'none')
cMap_cluster = pink(size(cellCountCluster_Area_prop, 1));
lumFactor = linspace(0.3, 1.8, 10);
for iB = 1:10
bb(iB).FaceColor = 'flat';
bb(iB).CData = cMap_cluster(iB,:); 
% bb(iB).CData = cMap_Area_MLfirst*lumFactor(iB);
end

% add lines across bars
cumYData = cumsum(cat(1, bb.YData));
for iB = 1:10
    xC = [bb(iB).XData(1)+bb(iB).BarWidth/2, bb(iB).XData(2)-bb(iB).BarWidth/2, bb(iB).XData(2)+bb(iB).BarWidth/2, ...
        bb(iB).XData(3)-bb(iB).BarWidth/2, bb(iB).XData(3)+bb(iB).BarWidth/2, bb(iB).XData(4)-bb(iB).BarWidth/2];
    yC = cumYData(iB, [1 2 2 3 3 4]);
    hold on
    line(reshape(yC, 2, 3), reshape(xC, 2, 3), 'Color', 'k', 'LineWidth', 0.5)
end

line([0.9 1], [5 5], 'Color', 'k', 'LineWidth', 3)
set(gca, 'YColor', 'none')
% print(gcf, fullfile(dirFig, 'cellGroupEachArea_barh_connected_nolabel'), '-r200', '-dtiff')
% print(gcf, fullfile(dirFig, 'cellGroupEachArea_barh_connected_nolabel'), '-depsc')

% color legend
fig_colorCluster = figure;
set(fig_colorCluster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 300 20 200]);
ax = subplot('Position', [0 0 1 1]);
image([1:10]');
colormap(pink(10));
axis off
line(repmat(get(gca, 'XLim')', 1, 11), repmat(0.5:1:10.5, 2, 1), 'Color', 'k')
set(gca, 'XTick', [], 'YTick', [])
box on
% print(fig_colorCluster, fullfile(dirFig, 'colormap_pink10_cluster'), '-depsc')
% print(fig_colorCluster, fullfile(dirFig, 'colormap_pink10_cluster'), '-r200', '-dtiff')


