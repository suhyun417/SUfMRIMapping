% playtimebuddy.m

clear all;

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat') %Clustering_CorrMap_4FPs_Movie123_ArtRHROI_probability.mat')


%%
setK = paramClustering_global.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS(:,iK) = sum(Clustering_meanROI.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
end

totalSS = Clustering_meanROI.totalSS_SU;
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;


%% Fig 3D: 2-D MDS plot showing K-means clustering results
D = pdist(Clustering_meanROI.matR, 'euclidean');
[Y2,stress,disparities] = mdscale(D,2);
[Y3,stress,disparities] = mdscale(D,3);

curK = 9; %11; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_meanROI.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

numROI = length(Clustering_meanROI.nameROI);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 200 910 675])

sp1 = subplot('Position', [0.15 0.2 0.8 0.7]);
imagesc(Clustering_meanROI.matR(indSortChan, :)')
set(sp1, 'CLim', [-1 1].*0.7)
set(sp1, 'YTick', 1:numROI, 'YTickLabel', Clustering_meanROI.nameROI)
locDiff = cat(1, find(diff(sortedClust)>0), length(sortedClust));
set(sp1, 'XTick', locDiff)
title(sprintf('Clustered cells from 4 FPs using max corr for each ROI: K=%d', curK))
xlabel('Cumulative number of cells')
colorbar;

sp2 = subplot('Position', [0.15 0.05 0.8 0.05]);
imagesc(Clustering_meanROI.catAreaID(indSortChan)')
set(sp2, 'YTick', 1, 'YTickLabel', 'Area info for each cell')
set(sp2, 'XTick', locDiff, 'XTickLabel', cat(1, locDiff(1), diff(locDiff)))
xlabel('Number of cells in each cluster')
colorbar;

% matIndClust_SU = cat(2, Clustering_moviemask_valid.resultKMeans.SU_indCluster); %cat(2, Clustering.resultKMeans.SU_indCluster);
% curK = 7; %
% [sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
indNewCluster = [4 2 3 6 7 5 1]; % [4 2 5 1 7 6 3]; %[1 2 3 4 5]; %[1 3 6 4 2 5]; %[2 4 1 3]; %[1 2 3 4 5]; %[3 4 2 5 1]; %[1 2 3 4 5 6 7]; %[4 1 6 3 5 2 7]; % cluster #4 is now cluster 1

% 

catChanID = cat(1, paramClustering.validChanID);
indMonkey = str2num(catChanID(:,5));

fig3a2=figure;
set(fig3a2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 700 700])
for iC=1:curK
    
    iK = indNewCluster(iC);
    
    curIndCells = resultProbClustering(iK).validIndCells;
    combCells = nchoosek(curIndCells,2);

    figure(fig3a2);
    line([Y2(combCells(:,1), 1) Y2(combCells(:,2),1)]', [Y2(combCells(:,1), 2) Y2(combCells(:,2), 2)]', 'Color', cMap(iC,:), 'LineWidth', 2);
    hold on;
    for iCell = 1:length(curIndCells)
        plot(Y2(curIndCells(iCell), 1), Y2(curIndCells(iCell), 2), 'o-', 'Marker', marker{indMonkey(curIndCells(iCell))},...
            'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(iC,:), 'LineWidth', 2)
        text(Y2(curIndCells(iCell),1)+1, Y2(curIndCells(iCell),2), catChanID(curIndCells(iCell),:))
        
%         plot3(Y3(curIndCells(iCell),1), Y3(curIndCells(iCell),2), Y3(curIndCells(iCell),3), 'o-','Marker', marker{indMonkey(curIndCells(iCell))}, ...
%             'LineWidth', 2, 'MarkerSize', 12,...
%             'Color', cMap(iC,:), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(iC,:));
%         text(Y3(curIndCells(iCell),1)+1, Y3(curIndCells(iCell),2), Y3(curIndCells(iCell),3), catChanID(curIndCells(iCell),:))
        hold on
    end    
%     figure(fig3a2); 
% %     curChan = indSortChan(sortedClust==iK);
%     p = plot(Y2(curIndCells, 1), Y2(curIndCells, 2), 'o','LineWidth', 2, 'MarkerSize', 10,...
%         'MarkerEdgeColor','k', 'MarkerFaceColor',cMap(iC,:), 'Color', cMap(iC,:)); %,...
%         %'Marker', marker{indMonkey(resultProbClustering(iC).validIndCells)});
%     hold on;
% %     plot(Y2(resultProbClustering(iC).validIndCells,1), Y2(resultProbClustering(iC).validIndCells,2), 'o','LineWidth', 2, 'MarkerSize', 10,...
% %         'MarkerEdgeColor','k', 'MarkerFaceColor', cMap(iC,:), 'Color', cMap(iC,:));
% %     text(Y2(curChan,1)+1, Y2(curChan,2), catChanID(curChan,:)) %paramCorr.validChanID(curChan,:))
%     hold on;
end

axis square

xlim([-32 35]) %xlim([-35 35]) %xlim([-25 30])
ylim([-10 11]) %ylim([-12 12])
set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 2,  'Box', 'off')
% print(fig3a2, fullfile(dirFig, 'figS4_MDS_6cluster'), '-depsc')