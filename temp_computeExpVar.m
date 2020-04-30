
curClustering = ClusteringSDF_TR; % ClusteringSDF_MION; % ClusteringSDF; %ClusteringSDFnorm;

setK = curClustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    indClust = curClustering.resultKMeans(iK).SU_indCluster; %Clustering_moviemask.resultKMeans(iK).SU_indCluster; %Clustering.resultKMeans(iK).SU_indCluster;
    [sortedClust, indSortedChan]=sort(indClust);
    
%     tExpVar=[];
%     for ii = 1:curK
%         tExpVar(ii,1) = Clustering.resultKMeans(iK).SU_sumD(ii)/(2*sum(sortedClust==ii));
%     end
    
    matWSS(iK,1) = sum(curClustering.resultKMeans(iK).SU_sumD); %sum(Clustering_moviemask.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
%     matExpVar(iK,1) = sum(tExpVar);

end

totalSS = curClustering.totalSS;
% [a, c, totalSS] = kmeans(matR_SU_moviemask', 1); %kmeans(matR_SU', 1);
betweenSS = totalSS-matWSS;
% totalVar = totalSS/(2*size(matR_SU,2)); %totalD/(2*size(matR_SU,2));

propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;


fig3b=figure;
set(fig3b, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 340 340])

curK=5;
plot(setK, propExplained, 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8); hold on
plot(curK, propExplained(curK-1), 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 8)
xlim([2 12])
ylim([0.5 0.9])
set(gca, 'TickDir', 'out', 'LineWidth', 2, 'Box', 'off', 'TickLength', [.025 .05])
set(gca, 'YTick', 0.5:.1:.9)