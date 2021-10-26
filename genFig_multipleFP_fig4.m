% genFig_multipleFP_fig4.m
%
% 2020/09/11 SHP
% Multiple face patch partially overlapping population
%   making figure 4B: bar graph of average correlation for each ROI of average map for each recording site

clear all;

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_probability.mat', 'Clustering_brainmask', 'param*') 
load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat', 'Clustering_meanROI')


%% Color maps
% ROIs
fid = fopen('/procdata/parksh/_macaque/Art/ROIs/i64_colorscale.pal');
A = fscanf(fid, '%s');
fclose(fid);
matRGB = sscanf(A(8:end), '#%2x%2x%2x', [3 inf])';
matRGB = flipud(unique(matRGB, 'rows', 'stable'));
matRGB = matRGB./255;

orderROI = [1 2 22 3 4 35 34 12 13 14 29 30 6 7 8 36 9 10 11 32 15 5 23 26 37 27 28 16 17 33 18 19 20 21 31 24 25]; % 1:37;

% Recording sites (4 face patches)
cMap_Area = [91 148 203; 237 28 35; 248 148 29; 6 177 102]./255; % from Kenji's schematic
cMap_Area(4, :) = cMap_Area(4, :).*0.7; % make the green a bit darker
orderArea = [4 1 2 3]; %ML-AF-AM-AAM

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


%% Fig 4A: composition of cell groups in each recording site
% Get the clustering results from voxel-based clustering
curK = 10; %13; % 9; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

cMap_cluster = pink(curK);

cellCountCluster_Area = NaN(curK, length(Clustering_brainmask.infoCells.setArea));
for iArea = 1:length(Clustering_brainmask.infoCells.setArea)
    compArea = [];
    compArea = sortedClust(ismember(indSortChan, find(Clustering_brainmask.infoCells.catAreaID == iArea)));
    cellCountCluster_Area(:, iArea) = histc(compArea, 1:curK);
end
cellCountCluster_Area_prop = cellCountCluster_Area./repmat(sum(cellCountCluster_Area), curK, 1);

reorderCluster = [8 1 9 5 3 7 4 10 6 2]; % ORder cell groups based on the number of neurons in each group (from largest to smallest)

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 200 1000 290]);
bb = barh(cellCountCluster_Area_prop(reorderCluster, orderArea)', 'stacked'); %barh(cellCountCluster_Area_prop(:, orderArea)', 'stacked');
set(gca, 'TickDir', 'out', 'Box', 'off');
set(gca, 'YDir', 'reverse')
set(gca, 'YTickLabel', Clustering_brainmask.infoCells.setArea(orderArea))
set(gca, 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XColor', 'none')

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
% print(gcf, fullfile(dirFig, 'cellGroupEachArea_barh_connected_nolabel_GroupReordered'), '-r200', '-dtiff')
% print(gcf, fullfile(dirFig, 'cellGroupEachArea_barh_connected_nolabel_GroupReordered'), '-depsc')

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

%% Principal components
aa = bsxfun(@minus, Clustering_meanROI.matR, mean(Clustering_meanROI.matR));
[coeff, score, latent] = pca(aa);

fig4B = figure;
set(fig4B, 'Color', 'w', 'PaperPositionMode', 'auto');
for iB = 1:3
sp(iB) = subplot(1,3,iB);
b(iB) = barh(coeff(orderROI, iB));
b(iB).FaceColor = 'flat';
    b(iB).CData = matRGB(orderROI, :);
    b(iB).BarWidth = 1;
    b(iB).EdgeColor = 'none';
    title(sp(iB), sprintf('PC %d (%2.2f%%)', iB, latent(iB)*100))
end
set(sp(:), 'XLim', [-0.4 0.4])
set(sp(:), 'YDir', 'reverse')
set(sp(1), 'YTick', 1:37, 'YTickLabel', Clustering_meanROI.nameROI(orderROI))
set(sp(2:3), 'YTick', 1:37, 'YTickLabel', [])
set(sp(:), 'TickDir', 'out', 'Box', 'off')

% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% set(b(:), 'FaceColor', 'flat', 'EdgeColor', 'none', 'CData', matRGB)


%
fig_scatter = figure;
set(fig_scatter, 'Color', 'w', 'PaperPositionMode', 'auto');
spp(1) = subplot(1,3,1);
scatter(score(:,1), score(:,2), 30, cMap_Area(Clustering_meanROI.catAreaID, :), 'fill');

spp(2) = subplot(1,3,2);
for iK = 1:curK
    scatter(score(sortedClust==iK,1), score(sortedClust==iK,2), 30, repmat(cMap_cluster(iK, :), sum(sortedClust==iK), 1), 'fill');
    hold on;
end

    


    
    
    % for iArea = 1:length(Clustering_meanROI.setArea)
%     matR_ROI_area{iArea} = mean(Clustering_meanROI.matR(Clustering_meanROI.catAreaID==iArea, orderROI));
% end
% 
% setArea_order = [4 1 2 3]; % ML-AF-AM-AAM
% 
% fig4A = figure;
% set(fig4A, 'Color', 'w', 'PaperPositionMode', 'auto');
% for iArea = 1:4
%     idArea = setArea_order(iArea);
%     sp(iArea) = subplot(4, 1, iArea);
%     b = bar(matR_ROI_area{idArea});
%     b.FaceColor = 'flat';
%     b.CData = matRGB(orderROI, :);
%     b.BarWidth = 1;
%     b.EdgeColor = 'none';
% end
% set(sp(:), 'Box', 'off', 'TickDir', 'out')
% set(sp(:), 'XColor', 'none')
% axis(sp(:), 'tight')

