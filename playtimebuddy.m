% playtimebuddy.m

clear all;

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

% load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability_pcares.mat') 
load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat')


[coeff, score, latent, tsquared, explained] = pca(zscore(Clustering_meanROI.matR));
figure
x = 1:37;
barh(x, coeff(:,1), 0.3, 'FaceColor', 'b', 'EdgeColor', 'none')
hold on
barh(x+0.2, coeff(:,2), 0.3, 'FaceColor', [0 0.7 0.7], 'EdgeColor', [0 0.7 0.7])
hold on
barh(x+0.4, coeff(:,3), 0.3, 'FaceColor', [0.7 0 0.7], 'EdgeColor', [0.7 0 0.7])
set(gca, 'YTick', 1:37, 'YTickLabel', Clustering_meanROI.nameROI)
xlabel('Coefficient from PCA')
hold off
legend(sprintf('PC1 (%d%%)', round(explained(1))), sprintf('PC2 (%d%%)', round(explained(2))), sprintf('PC3 (%d%%)', round(explained(3))))
title('Before remove first PC from fMRI')

%%
setK = paramClustering_global.setK; %Clustering.setK;

matWSS=[];matWSS_roi=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS(:,iK) = sum(Clustering_meanROI.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
    matWSS_roi(:,iK) = sum(Clustering_meanROI.resultKMeans(iK).roi_sumD); %
end

totalSS = Clustering_meanROI.totalSS_SU;
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(setK, propExplained'.*100, 'ko-', 'MarkerFaceColor', 'w'); hold on
xlabel('Number of cluster (K)')
ylabel('Explained variance (%)')
title('Clustering using mean r for each ROI')
set(gca, 'XTick', setK)
set(gca, 'TickDir', 'out', 'Box', 'off')

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(setK(1:end-1), diff(propExplained').*100)
hold on
plot(setK(1:end-1), mean(diff(propExplained').*100, 2), 'ko-', 'LineWidth', 2)
title('difference of explained variance for each K: using mean ROI')

%% Colormap
fname = '/procdata/parksh/_macaque/Art/Anatomy/_suma/BCWYRColorMap.txt';
ttt = dlmread(fname);
cMap_corrSUMA = ttt(:, 1:3);
clear ttt


%% Fig 3D: 2-D MDS plot showing K-means clustering results
% D = pdist(Clustering_meanROI.matR, 'euclidean');
% [Y2,stress,disparities] = mdscale(D,2);
% [Y3,stress,disparities] = mdscale(D,3);

curK = 9; %13; % 9; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_meanROI.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

numROI = length(Clustering_meanROI.nameROI);
orderROI = [1 2 22 3 4 35 34 12 13 14 29 30 6 7 8 36 9 10 11 32 15 5 23 26 37 27 28 16 17 33 18 19 20 21 31 24 25]; % 1:37;
cMap_Area = [91 148 203; 237 28 35; 248 148 29; 6 177 102]./255; % from Kenji's schematic

%
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 200 910 675])

% % make blue-white-red colorbar
% cval = 0.5;
% cmin = -cval; cmax = cval;
% colornum = 256;
% colorInput = [1 0 0; 1 1 1; 0 0 1];
% oldSteps = linspace(-1, 1, length(colorInput));
% newSteps = linspace(-1, 1, colornum);
% for j=1:3 % RGB
%     newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
% end
% endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
% newmap = squeeze(newmap_all(1:endPoint, :));
% % figure(gcf)
% % set(gca, 'CLim', [cmin cmax])
% % colormap(flipud(newmap))

sp1 = subplot('Position', [0.15 0.2 0.8 0.7]);
imagesc(Clustering_meanROI.matR(indSortChan, orderROI)')
set(sp1, 'CLim', [-1 1].*0.5)
set(sp1, 'YTick', 1:numROI, 'YTickLabel', Clustering_meanROI.nameROI(orderROI))
locDiff = cat(1, find(diff(sortedClust)>0), length(sortedClust));
set(sp1, 'XTick', locDiff+0.5, 'XTickLabel', locDiff)
title(sprintf('Clustered cells from 4 FPs using mean corr for each ROI: K=%d', curK))
xlabel('Cumulative number of cells')
line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k')
colormap(sp1, cMap_corrSUMA)
set(gca, 'TickDir', 'out', 'Box', 'off')
colorbar;

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 200 910 70])
sp2 = subplot('Position', [0.15 0.05 0.8 0.8]);
imagesc(Clustering_meanROI.catAreaID(indSortChan)')
set(sp2, 'YTick', 1, 'YTickLabel', 'Area info for each cell')
set(sp2, 'XTick', locDiff+0.5, 'XTickLabel', cat(1, locDiff(1), diff(locDiff)))
line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k')
xlabel('Number of cells in each cluster')
colormap(sp2, cMap_Area)
colorbar;
axis off



%
cellCountCluster_Area = NaN(curK, length(Clustering_meanROI.setArea));
for iSubj = 1:length(Clustering_meanROI.setArea)
    compArea = [];
    compArea = sortedClust(ismember(indSortChan, find(Clustering_meanROI.catAreaID == iSubj)));
    cellCountCluster_Area(:, iSubj) = histc(compArea, 1:curK);
end
cellCountCluster_Area_prop = cellCountCluster_Area./repmat(sum(cellCountCluster_Area), curK, 1);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 700 1720 360]);
for iSubj = 1:length(Clustering_meanROI.setArea)
    sp(iSubj) = subplot(1, length(Clustering_meanROI.setArea), iSubj);
    pie(sp(iSubj), cellCountCluster_Area(:, iSubj));
    title(Clustering_meanROI.setArea{iSubj})
end

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 700 355 730]);
for iSubj = 1:length(Clustering_meanROI.setArea)
    sp(iSubj) = subplot(length(Clustering_meanROI.setArea), 1, iSubj);
    bar(sp(iSubj), cellCountCluster_Area_prop(:, iSubj).*100);
    title(Clustering_meanROI.setArea{iSubj})
end
xlabel(sp(4), 'Cluster ID')
ylabel(sp(2), 'Percent of cells in each cluster (%)')

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
print(fig_colorROI, fullfile(dirFig, 'colormap_ROIset01_reorder'), '-depsc')
print(fig_colorROI, fullfile(dirFig, 'colormap_ROIset01_reorder'), '-r200', '-dtiff')

%% For max r for each ROI
setK = paramClustering_global.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS(:,iK) = sum(Clustering_maxabsROI.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
end

totalSS = Clustering_maxabsROI.totalSS_SU;
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(setK, propExplained'.*100, 'ko-'); hold on
xlabel('Number of cluster (K)')
ylabel('Explained variance (%)')
title('Clustering using max r for each ROI')
set(gca, 'XTick', setK)

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(setK(1:end-1), diff(propExplained').*100)
hold on
plot(setK(1:end-1), mean(diff(propExplained').*100, 2), 'ko-', 'LineWidth', 2)
title('difference of explained variance for each K: using max ROI')

%% Fig 3D: 2-D MDS plot showing K-means clustering results
D = pdist(Clustering_maxabsROI.matR, 'euclidean');
[Y2,stress,disparities] = mdscale(D,2);
[Y3,stress,disparities] = mdscale(D,3);

curK = 5; %8; %11; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_maxabsROI.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

numROI = length(Clustering_maxabsROI.nameROI);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 200 910 675])

sp1 = subplot('Position', [0.15 0.2 0.8 0.7]);
imagesc(Clustering_maxabsROI.matR(indSortChan, :)')
set(sp1, 'CLim', [-1 1].*0.7)
set(sp1, 'YTick', 1:numROI, 'YTickLabel', Clustering_maxabsROI.nameROI)
locDiff = cat(1, find(diff(sortedClust)>0), length(sortedClust));
set(sp1, 'XTick', locDiff)
title(sprintf('Clustered cells from 4 FPs using max corr for each ROI: K=%d', curK))
xlabel('Cumulative number of cells')
colorbar;

sp2 = subplot('Position', [0.15 0.05 0.8 0.05]);
imagesc(Clustering_maxabsROI.catAreaID(indSortChan)')
set(sp2, 'YTick', 1, 'YTickLabel', 'Area info for each cell')
set(sp2, 'XTick', locDiff, 'XTickLabel', cat(1, locDiff(1), diff(locDiff)))
xlabel('Number of cells in each cluster')
colorbar;


%
cellCountCluster_Area = NaN(curK, length(Clustering_maxabsROI.setArea));
for iSubj = 1:length(Clustering_maxabsROI.setArea)
    compArea = [];
    compArea = sortedClust(ismember(indSortChan, find(Clustering_maxabsROI.catAreaID == iSubj)));
    cellCountCluster_Area(:, iSubj) = histc(compArea, 1:curK);
end
cellCountCluster_Area_prop = cellCountCluster_Area./repmat(sum(cellCountCluster_Area), curK, 1);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 700 1720 360]);
for iSubj = 1:length(Clustering_maxabsROI.setArea)
    sp(iSubj) = subplot(1, length(Clustering_maxabsROI.setArea), iSubj);
    pie(sp(iSubj), cellCountCluster_Area(:, iSubj));
    title(Clustering_maxabsROI.setArea{iSubj})
end

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 700 355 730]);
for iSubj = 1:length(Clustering_maxabsROI.setArea)
    sp(iSubj) = subplot(length(Clustering_maxabsROI.setArea), 1, iSubj);
    bar(sp(iSubj), cellCountCluster_Area_prop(:, iSubj).*100);
    title(Clustering_maxabsROI.setArea{iSubj})
end
xlabel(sp(4), 'Cluster ID')
ylabel(sp(2), 'Percent of cells in each cluster (%)')



%% Clustering results using all the voxels, not ROI
clear all;

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat', 'Clustering_meanROI')
load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_probability.mat', 'Clustering_brainmask', 'param*') 

% [coeff, score, latent, tsquared, explained] = pca(zscore(Clustering_brainmask.matR));
% figure
% x = 1:37;
% barh(x, coeff(:,1), 0.3, 'FaceColor', 'b', 'EdgeColor', 'none')
% hold on
% barh(x+0.2, coeff(:,2), 0.3, 'FaceColor', [0 0.7 0.7], 'EdgeColor', [0 0.7 0.7])
% hold on
% barh(x+0.4, coeff(:,3), 0.3, 'FaceColor', [0.7 0 0.7], 'EdgeColor', [0.7 0 0.7])
% set(gca, 'YTick', 1:37, 'YTickLabel', Clustering_meanROI.nameROI)
% xlabel('Coefficient from PCA')
% hold off
% legend(sprintf('PC1 (%d%%)', round(explained(1))), sprintf('PC2 (%d%%)', round(explained(2))), sprintf('PC3 (%d%%)', round(explained(3))))
% title('Before remove first PC from fMRI')

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


%% Correlation matrix

curK = 10; %13; % 9; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

% Colormap
cMap_Area = [91 148 203; 237 28 35; 248 148 29; 6 177 102]./255; % from Kenji's schematic
cMap_Area(4, :) = cMap_Area(4, :).*0.7;
orderArea = [4 1 2 3]; %ML-AF-AM-AAM
cMap_Area_MLfirst = cMap_Area(orderArea, :);

fname = '/procdata/parksh/_macaque/Art/Anatomy/_suma/BCWYRColorMap.txt';
ttt = dlmread(fname);
cMap_corrSUMA = ttt(:, 1:3);
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
%
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 200 910 675])

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

% colorbar;
set(sp1, 'XTickLabel', [])
spp1 = axes('Position', sp1.Position, 'Color', 'none', 'XAxisLocation', 'top');
spp1.XLim = sp1.XLim;
set(spp1, 'XTick', locDiff+0.5, 'XTickLabel', locDiff, 'TickDir', 'out', 'XTickLabel', [], 'Box', 'off', 'YColor', 'none')

sp2 = subplot('Position', [0.2 0.08 0.7 0.05]);
imagesc(Clustering_brainmask.infoCells.catAreaID(indSortChan_reorder)') %(indSortChan)')
set(sp2, 'YTick', []) %1, 'YTickLabel', 'Area info for each cell')
set(sp2, 'XTick', locDiff+0.5, 'XTickLabel', [], 'TickDir', 'out'); %cat(1, locDiff(1), diff(locDiff)))
line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k')
% xlabel('Number of cells in each cluster')
colormap(sp2, cMap_Area)
spp2 = axes('Position', sp2.Position, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
spp2.Clipping = 'off';
spp2.XLim = sp2.XLim;
spp2.YLim = sp2.YLim;
line([locDiff+0.5 locDiff+0.5]', repmat([0.5;3], 1, length(locDiff)), 'Color', 'k')
% colorbar;


% sp1 = subplot('Position', [0.15 0.2 0.8 0.7]);
% imagesc(Clustering_brainmask.matR(:, indSortChan))
% set(sp1, 'CLim', [-1 1].*0.5)
% % set(sp1, 'YTick', 1:numROI, 'YTickLabel', Clustering_brainmask.nameROI(orderROI))
% locDiff = cat(1, find(diff(sortedClust)>0), length(sortedClust));
% set(sp1, 'XTick', locDiff)
% title(sprintf('Clustered cells from 4 FPs using voxels within brain: K=%d', curK))
% xlabel('Cumulative number of cells')
% colorbar;
% 
% % figure;
% % set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 200 910 675])
% sp2 = subplot('Position', [0.15 0.05 0.8 0.05]);
% imagesc(Clustering_brainmask.infoCells.catAreaID(indSortChan)')
% set(sp2, 'YTick', 1, 'YTickLabel', 'Area info for each cell')
% set(sp2, 'XTick', locDiff, 'XTickLabel', cat(1, locDiff(1), diff(locDiff)))
% xlabel('Number of cells in each cluster')
% % colormap(sp2, cMap_Area)
% colorbar;


%
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
print(gcf, fullfile(dirFig, 'cellGroupEachArea_barh_connected_nolabel'), '-r200', '-dtiff')
print(gcf, fullfile(dirFig, 'cellGroupEachArea_barh_connected_nolabel'), '-depsc')
% print(gcf, fullfile(dirFig, 'ClusterComposition_KMeansClustering_brainmaskVoxels_K10_forEachArea'), '-depsc')

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
print(fig_colorCluster, fullfile(dirFig, 'colormap_pink10_cluster'), '-depsc')
print(fig_colorCluster, fullfile(dirFig, 'colormap_pink10_cluster'), '-r200', '-dtiff')



figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 700 1720 360]);
for iSubj = 1:length(Clustering_brainmask.infoCells.setArea)
    sp(iSubj) = subplot(1, length(Clustering_brainmask.infoCells.setArea), iSubj);
    pie(sp(iSubj), cellCountCluster_Area(:, iSubj));
    title(Clustering_brainmask.infoCells.setArea{iSubj})
end

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [300 700 355 730]);
for iSubj = 1:length(Clustering_brainmask.infoCells.setArea)
    sp(iSubj) = subplot(length(Clustering_brainmask.infoCells.setArea), 1, iSubj);
    bar(sp(iSubj), cellCountCluster_Area_prop(:, iSubj).*100);
    title(Clustering_brainmask.infoCells.setArea{iSubj})
end
xlabel(sp(4), 'Cluster ID')
ylabel(sp(2), 'Percent of cells in each cluster (%)')


%% mean corrleation per ROI for Figure 2 example cells
setExampleCellIDs = {'33Dav', '130AFMoc', '097aMat', '10Dan', '27Dav', '065aTor', '39AMWas', '117AMMoc', '25Dav', '022bSpi', '51AMWas', '109AMMoc', '16Dav', '045aSpi', '33AMWas', '05Dan'};

figure;
for iCell = 1:length(setExampleCellIDs)
subplot(4, 4, iCell);
curCellID = setExampleCellIDs{iCell};
bar(Clustering_meanROI.matR(ismember(Clustering_meanROI.catChanID, curCellID), orderROI))
end

fid = fopen('/procdata/parksh/_macaque/Art/ROIs/i64_colorscale.pal');
A = fscanf(fid, '%s');
fclose(fid);
matRGB = sscanf(A(8:end), '#%2x%2x%2x', [3 inf])';
matRGB = flipud(unique(matRGB, 'rows', 'stable'));
matRGB = matRGB./255;

% 2018a bar color
figure;
for iCell = 1:length(setExampleCellIDs)
sp(iCell) = subplot(4, 4, iCell);
curCellID = setExampleCellIDs{iCell};
b = bar(Clustering_meanROI.matR(ismember(Clustering_meanROI.catChanID, curCellID), orderROI));
b.FaceColor = 'flat';
b.CData = matRGB(orderROI, :);
b.EdgeColor = 'none';
end
set(sp(1:4), 'YLim', [-0.05 0.4])
set(sp(5:8), 'YLim', [-0.5 0.1])
set(sp(9:12), 'YLim', [-0.1 0.4])
set(sp(13:16), 'YLim', [-0.2 0.4])
set(sp(:), 'TickDir', 'out', 'Box', 'off', 'XTick', [])
set(gcf, 'COlor', 'w')
set(sp(:), 'XColor', [1 1 1], 'YColor', [0 0 0])

% 2018a bar color
figBar = figure;
set(figBar, 'Color', 'w', 'PaperPositionMode', 'auto')
for iCell = 1:length(setExampleCellIDs)
sp(iCell) = subplot(4, 4, iCell);
curCellID = setExampleCellIDs{iCell};
b = bar(Clustering_meanROI.matR(ismember(Clustering_meanROI.catChanID, curCellID), orderROI), 'hist');
b.EdgeColor = 'none';
b.FaceAlpha = 0.6;
end
set(sp(1:4), 'YLim', [-0.05 0.4])
set(sp(5:8), 'YLim', [-0.5 0.1])
set(sp(9:12), 'YLim', [-0.1 0.4])
set(sp(13:16), 'YLim', [-0.2 0.4])
set(sp(:), 'TickDir', 'out', 'Box', 'off', 'XTick', [])
set(gcf, 'COlor', 'w')
set(sp(:), 'XColor', [1 1 1], 'YColor', [0 0 0])


%%
load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_probability.mat', 'Clustering_brainmask', 'param*')
D_corr_voxel = pdist(Clustering_brainmask.matR', 'spearman');
matD_corr_voxel = squareform(D_corr_voxel);
figure
imagesc(matD_corr_voxel)
matPairWiseCorr_corrmapspace = 1-matD_corr_voxel;

% within- & between-area comparison
for iArea = 1:length(Clustering_brainmask.infoCells.setArea)
    indCellArea{iArea} = Clustering_brainmask.infoCells.catAreaID==iArea;
end

for iArea1 = 1:length(Clustering_brainmask.infoCells.setArea)
    for iArea2 = 1:length(Clustering_brainmask.infoCells.setArea)
        
        resultsR_SUmap_area(iArea1, iArea2).nameArea = {Clustering_brainmask.infoCells.setArea{iArea1}, Clustering_brainmask.infoCells.setArea{iArea2}};
        resultsR_SUmap_area(iArea1, iArea2).matR = matPairWiseCorr_corrmapspace(indCellArea{iArea1}, indCellArea{iArea2});
        if iArea1==iArea2
            tempL = tril(resultsR_SUmap_area(iArea1, iArea2).matR, -1);
            resultsR_SUmap_area(iArea1, iArea2).vectR = tempL(abs(tempL)>0);
        else
            resultsR_SUmap_area(iArea1, iArea2).vectR = resultsR_SUmap_area(iArea1, iArea2).matR(:);
        end
        
        resultsR_SUmap_area(iArea1, iArea2).meanR = mean(resultsR_SUmap_area(iArea1, iArea2).vectR);
        resultsR_SUmap_area(iArea1, iArea2).medianR = median(resultsR_SUmap_area(iArea1, iArea2).vectR);
        resultsR_SUmap_area(iArea1, iArea2).steR = std(resultsR_SUmap_area(iArea1, iArea2).vectR)/length(resultsR_SUmap_area(iArea1, iArea2).vectR);
    end
end
         
aa = struct2cell(resultsR_SUmap_area);
mat_meanR = squeeze(cell2mat(aa(4, :, :)));


figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [200 200 850 600])
edges = -1:0.2:1;
for iArea1 = 1:4
    for iArea2 = 1:4
        spp(iArea1, iArea2) = subplot(4,4,4*(iArea1-1)+iArea2);
        h = histogram(resultsR_SUmap_area(iArea1, iArea2).vectR, edges);        
        h.Normalization = 'probability';
        h.EdgeColor = 'none';
        hold on
        line(zeros(1, 2), get(gca, 'YLim'), 'Color', 'k');
        line(ones(1, 2).*resultsR_SUmap_area(iArea1, iArea2).medianR, get(gca, 'YLim'), 'Color', 'r', 'LineWidth', 2);
    end
end
set(spp(:), 'TickDir', 'out', 'Box', 'off')
set(spp(:), 'XColor', [0 0 0], 'YColor', [0 0 0])




% within- & between-animal comparison
setSubj = unique(Clustering_brainmask.infoCells.catSubjID);
for iSubj = 1:length(setSubj)
    indCellSubj{iSubj} = Clustering_brainmask.infoCells.catSubjID==setSubj(iSubj);
end

for iSubj1 = 1:length(setSubj)
    for iSubj2 = 1:length(setSubj)
        
        resultsR_SUmap_subj(iSubj1, iSubj2).subjID = [setSubj(iSubj1), setSubj(iSubj2)];
        resultsR_SUmap_subj(iSubj1, iSubj2).matR = matPairWiseCorr_corrmapspace(indCellSubj{iSubj1}, indCellSubj{iSubj2});
        if iSubj1==iSubj2
            tempL = tril(resultsR_SUmap_subj(iSubj1, iSubj2).matR, -1);
            resultsR_SUmap_subj(iSubj1, iSubj2).vectR = tempL(abs(tempL)>0);
        else
            resultsR_SUmap_subj(iSubj1, iSubj2).vectR = resultsR_SUmap_subj(iSubj1, iSubj2).matR(:);
        end
        
        resultsR_SUmap_subj(iSubj1, iSubj2).meanR = mean(resultsR_SUmap_subj(iSubj1, iSubj2).vectR);
        resultsR_SUmap_subj(iSubj1, iSubj2).medianR = median(resultsR_SUmap_subj(iSubj1, iSubj2).vectR);
        resultsR_SUmap_subj(iSubj1, iSubj2).steR = std(resultsR_SUmap_subj(iSubj1, iSubj2).vectR)/length(resultsR_SUmap_subj(iSubj1, iSubj2).vectR);
    end
end
         
tempCellSubj = struct2cell(resultsR_SUmap_subj);
mat_meanR_sub = squeeze(cell2mat(tempCellSubj(4, :, :)));

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [200 200 850 600])
edges = -1:0.2:1;
for iSubj1 = 1:10
    for iSubj2 = 1:10
        spp(iSubj1, iSubj2) = subplot(10,10,10*(iSubj1-1)+iSubj2);
        h = histogram(resultsR_SUmap_subj(iSubj1, iSubj2).vectR, edges);        
        h.Normalization = 'probability';
        h.EdgeColor = 'none';
        hold on
        line(zeros(1, 2), get(gca, 'YLim'), 'Color', 'k');
        line(ones(1, 2).*resultsR_SUmap_subj(iSubj1, iSubj2).medianR, get(gca, 'YLim'), 'Color', 'r', 'LineWidth', 2);
    end
end
set(spp(:), 'TickDir', 'out', 'Box', 'off')
set(spp(:), 'XColor', [0 0 0], 'YColor', [0 0 0])


%% making a map of areas showing high correlation with at least 10% of cells from all face patches
load('/procdata/parksh/_macaque/CorrMap_SU_AllCellsArt_corticalFPMerged.mat', 'corrMap_Area')
 
critCorr = 0.5; %0.5; %0.4; %35;
 
for iArea = 1:4
setR = corrMap_Area(iArea).matR;
matValidVox = abs(setR)>critCorr;
tMatVox{iArea} = sum(matValidVox, 2)./size(matValidVox,2);
end


nx = 40; ny = 64; nz = 32; %LR-AP-DV
critRatio = 0.01; %0.01; %0.03; %0.01; %0.1;
matVoxHighCorrCell = sum(cat(2, tMatVox{:})>critRatio, 2);
matVoxHighCorrCell_3D = reshape(matVoxHighCorrCell, [nx, ny, nz]);

ttt_s = permute(matVoxHighCorrCell_3D, [3 1 2]); %[3 2 1]); % %[3 1 2]: coronal sections %[2 1 3]: axial sections %[3 2 1]: sagittal sections
% fig_sections = figure;
for iS = 11:50
subplot(5, 8, iS-10);
imagesc(ttt_s(:,:,iS))
set(gca, 'CLim', [0 4])
colormap(parula(5))
title(sprintf('section %d', iS))
axis off
end

ttt_s = permute(matVoxHighCorrCell_3D, [2 1 3]); %[3 2 1]); % %[3 1 2]: coronal sections %[2 1 3]: axial sections %[3 2 1]: sagittal sections
fig_sections = figure;
for iS = 1:size(ttt_s, 3)
figure(fig_sections); clf
imagesc(ttt_s(:,:,iS))
set(gca, 'CLim', [0 4])
colormap(parula(5))
title(sprintf('section %d', iS))
colorbar;
input('')
end


ttt_s = permute(matVoxHighCorrCell_3D, [3 2 1]); % %[3 1 2]: coronal sections %[2 1 3]: axial sections %[3 2 1]: sagittal sections
fig_sections = figure;
for iS = 1:size(ttt_s, 3)
subplot(8, 5, iS);
imagesc(ttt_s(:,:,iS))
set(gca, 'CLim', [0 4])
colormap(parula(5))
title(sprintf('section %d', iS))
axis off
end

numCol = 8;
numRow = 5;
marginCol = 0.01;
marginRow = 0.01;
widthSP = (1-marginCol*2)/numCol;
heightSP = (1-marginRow*2)/numRow;

% %%
% % matIndClust_SU = cat(2, Clustering_moviemask_valid.resultKMeans.SU_indCluster); %cat(2, Clustering.resultKMeans.SU_indCluster);
% % curK = 7; %
% % [sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
% indNewCluster = [4 2 3 6 7 5 1]; % [4 2 5 1 7 6 3]; %[1 2 3 4 5]; %[1 3 6 4 2 5]; %[2 4 1 3]; %[1 2 3 4 5]; %[3 4 2 5 1]; %[1 2 3 4 5 6 7]; %[4 1 6 3 5 2 7]; % cluster #4 is now cluster 1
% 
% % 
% 
% catChanID = cat(1, paramClustering.validChanID);
% indMonkey = str2num(catChanID(:,5));
% 
% fig3a2=figure;
% set(fig3a2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 700 700])
% for iC=1:curK
%     
%     iK = indNewCluster(iC);
%     
%     curIndCells = resultProbClustering(iK).validIndCells;
%     combCells = nchoosek(curIndCells,2);
% 
%     figure(fig3a2);
%     line([Y2(combCells(:,1), 1) Y2(combCells(:,2),1)]', [Y2(combCells(:,1), 2) Y2(combCells(:,2), 2)]', 'Color', cMap(iC,:), 'LineWidth', 2);
%     hold on;
%     for iCell = 1:length(curIndCells)
%         plot(Y2(curIndCells(iCell), 1), Y2(curIndCells(iCell), 2), 'o-', 'Marker', marker{indMonkey(curIndCells(iCell))},...
%             'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(iC,:), 'LineWidth', 2)
%         text(Y2(curIndCells(iCell),1)+1, Y2(curIndCells(iCell),2), catChanID(curIndCells(iCell),:))
%         
% %         plot3(Y3(curIndCells(iCell),1), Y3(curIndCells(iCell),2), Y3(curIndCells(iCell),3), 'o-','Marker', marker{indMonkey(curIndCells(iCell))}, ...
% %             'LineWidth', 2, 'MarkerSize', 12,...
% %             'Color', cMap(iC,:), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(iC,:));
% %         text(Y3(curIndCells(iCell),1)+1, Y3(curIndCells(iCell),2), Y3(curIndCells(iCell),3), catChanID(curIndCells(iCell),:))
%         hold on
%     end    
% %     figure(fig3a2); 
% % %     curChan = indSortChan(sortedClust==iK);
% %     p = plot(Y2(curIndCells, 1), Y2(curIndCells, 2), 'o','LineWidth', 2, 'MarkerSize', 10,...
% %         'MarkerEdgeColor','k', 'MarkerFaceColor',cMap(iC,:), 'Color', cMap(iC,:)); %,...
% %         %'Marker', marker{indMonkey(resultProbClustering(iC).validIndCells)});
% %     hold on;
% % %     plot(Y2(resultProbClustering(iC).validIndCells,1), Y2(resultProbClustering(iC).validIndCells,2), 'o','LineWidth', 2, 'MarkerSize', 10,...
% % %         'MarkerEdgeColor','k', 'MarkerFaceColor', cMap(iC,:), 'Color', cMap(iC,:));
% % %     text(Y2(curChan,1)+1, Y2(curChan,2), catChanID(curChan,:)) %paramCorr.validChanID(curChan,:))
% %     hold on;
% end
% 
% axis square
% 
% xlim([-32 35]) %xlim([-35 35]) %xlim([-25 30])
% ylim([-10 11]) %ylim([-12 12])
% set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 2,  'Box', 'off')
% % print(fig3a2, fullfile(dirFig, 'figS4_MDS_6cluster'), '-depsc')