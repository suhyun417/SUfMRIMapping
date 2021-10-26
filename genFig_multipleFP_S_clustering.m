% genFig_multipleFP_S_clustering.m
%
% 2021/04/29 SHP
% Supplemental figure for multiple face patch paper
% showing clustering related details (similar to the Fig. 3 of Park et al. 2017 Neuron)

%% Settings
ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS';
    dirProcdata = '/Volumes/PROCDATA';
    dirLibrary = '/Volumes/LIBRARY';
else % on virtual machine
    dirProjects = '/projects';
    dirProcdata = '/procdata';
    dirLibrary = '/library';
end
    
% Add necessary toolbox % Should be 
addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution
addpath(fullfile(dirProjects, 'parksh/_toolbox/afni_matlab'))
addpath(fullfile(dirProjects, 'parksh/_toolbox/hslcolormap'))

% Directory for saving figures as graphic files
dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

%% Load data
% load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat', 'Clustering_meanROI')
load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_probability.mat', 'Clustering_brainmask', 'param*') 



%% Explained variance elbow plot
setK = paramClustering_global.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS(:,iK) = sum(Clustering_brainmask.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
end

totalSS = Clustering_brainmask.totalSS_SU;
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

fig3b=figure;
set(fig3b, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 340 340])
curK=10;
% plot(setK, propExplained.*100, 'k-', 'LineWidth', 1); hold on
plot(setK, propExplained.*100, 'ko-', 'LineWidth', 1, 'MarkerFaceColor', 'w', 'MarkerSize', 8); hold on
plot(curK, propExplained(:, curK-1).*100, 'ko', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 8)
xlim([2 20])
ylim([35 75])
axis square
set(gca, 'TickDir', 'out', 'LineWidth', 2, 'Box', 'off', 'TickLength', [.025 .05], 'FontSize', 15)
set(gca, 'YTick', 35:10:75)

% save
print(fig3b, fullfile(dirFig, 'multipleFP_FigS_KMeansClustering_brainmaskVoxels_expVar'), '-depsc')

%% Difference & variance of explained variance across different Ks
deltaPropExplained = diff(mean(propExplained).*100);
fig3c = figure;
set(fig3c, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 340 340], 'Name', 'Difference of explained variance');
plot(setK(2:end), diff(mean(propExplained).*100), 'ko-', 'LineWidth', 1, 'MarkerFaceColor', 'w', 'MarkerSize', 8); hold on;
curK=10;
plot(curK, deltaPropExplained(curK-2), 'ko', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 8)
xlim([2 20])
axis square
set(gca, 'YTick', 0:3:12)
set(gca, 'TickDir', 'out', 'LineWidth', 2, 'Box', 'off', 'TickLength', [.025 .05], 'FontSize', 15)

% save
print(fig3c, fullfile(dirFig, 'multipleFP_FigS_KMeansClustering_brainmaskVoxels_expVar_delta'), '-depsc')

%
stdPropExplained = std(propExplained.*100);
fig3d = figure;
set(fig3d, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 340 340], 'Name', 'Variance of explained variance');
plot(setK, stdPropExplained, 'ko-', 'LineWidth', 1, 'MarkerFaceColor', 'w', 'MarkerSize', 8); hold on;
curK=10;
plot(curK, stdPropExplained(curK-1), 'ko', 'LineWidth', 1, 'MarkerFaceColor', 'k', 'MarkerSize', 8)
xlim([2 20])
ylim([0 0.18])
axis square
set(gca, 'YTick', 0:0.06:0.18)
set(gca, 'TickDir', 'out', 'LineWidth', 2, 'Box', 'off', 'TickLength', [.025 .05], 'FontSize', 15)

% save
print(fig3d, fullfile(dirFig, 'multipleFP_FigS_KMeansClustering_brainmaskVoxels_expVar_std'), '-depsc')


%% Fig 3C: Probability of the clustering
% K=10; %10; %7; %8; %6; %7;
% locMode = find(propExplained(:,K-1)==mode(propExplained(:,K-1)));
% locMin = find(propExplained(:,K-1)==min(propExplained(:,K-1)));
% indClust = Clustering_brainmask.resultKMeans(K-1).SU_indCluster(:, locMode(1)); %
% [sortClust, sortCell] = sort(Clustering_brainmask.resultKMeans(K-1).SU_indCluster(:, locMode(1)));
% 
% % orderROI = [4 2 5]; %[4 1 3 2 5]; %% 1: parafoveal:ventral stream / 2: face patches / 3: Eccentric visual cortex / 4: Central V1 / 5: MT
% orderClust_Cell = 1:K; %[1 2 5 6 3 7 4]; %1:K; % [1 2 5 6 3 7 4]; %[4 2 5 1 7 6 3]; %
% 
% reOrdSortCell=[];
% for iCC=1:K
%     reOrdSortCell = cat(1, reOrdSortCell, find(indClust==orderClust_Cell(iCC)));
% end
% 
% fig3c=figure;
% set(fig3c, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 510 470])
% imagesc(Clustering_brainmask.resultKMeans(K-1).matProb(reOrdSortCell, reOrdSortCell)); % how probable
% axis square
% axis off
% set(gca, 'CLim', [0 1])
% 
% % % save
% % print(fig3c, fullfile(dirFig, sprintf('figRev3C_Likelihood_%dMeans', K)), '-depsc')
% % print(fig3c, fullfile(dirFig, sprintf('figRev3C_Likelihood_%dMeans', K)), '-dtiff', '-r200')
% % % print(fig3c, fullfile(dirFig, sprintf('figRevS4_Likelihood_%dMeans', K)), '-dtiff', '-r200')


%% Projected to PCs
[coeff_e, score_e, latent_e, tsquared_e, explained_e] = pca(Clustering_brainmask.matR(:, locFaceCell)'); %, 'Economy', false);

curK=10; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

% 
cMap = pink(curK);
cMap_Area = [179 226 205; 141 160 203; 252 141 98; 231 41 138]./255; % 
catChanID = cat(1, Clustering_brainmask.infoCells.catChanID);
% indMonkey = str2num(catChanID(:,5));
% ordering of cells: ML-AF-AM-AAM, while maintaining the grouping
for iK = 1:curK
    tempS(iK).indSortChan_org = indSortChan(sortedClust==iK);
    tempS(iK).sortedClust_org = sortedClust(sortedClust==iK);
    tempS(iK).numCells = sum(sortedClust==iK);
end
reorderCluster = [8 1 9 5 3 7 4 10 6 2]; % ORder cell groups based on the number of neurons in each group (from largest to smallest)
indSortChan_reorder = cat(1, tempS(reorderCluster).indSortChan_org);
sortedClust_reorder = cat(1, tempS(reorderCluster).sortedClust_org);

% Colored by recording sites
fig3a2=figure;
set(fig3a2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 700 700])
Y2(:,1) = score_e(:,1); %.*(-1); % to match it to the MDS plot
Y2(:,2) = score_e(:,2);
for iC=1:curK
    
    iK = reorderCluster(iC);
    
    curIndCells = indSortChan_reorder(sortedClust_reorder==iK); %resultProbClustering(iK).validIndCells;
    combCells = nchoosek(curIndCells,2);

    figure(fig3a2);
%     line([Y2(combCells(:,1), 1) Y2(combCells(:,2),1)]', [Y2(combCells(:,1), 2) Y2(combCells(:,2), 2)]', 'Color', cMap(iC,:), 'LineWidth', 2);
    hold on;
    for iCell = 1:length(curIndCells)
        marker = 'o';
%         if ismember(curIndCells(iCell), locFaceCell)
%             marker = '^';
%         end
        plot(Y2(curIndCells(iCell), 1), Y2(curIndCells(iCell), 2), 'o-', 'Marker', marker, ... %marker{indMonkey(curIndCells(iCell))},...
            'MarkerSize', 10, 'MarkerEdgeColor', cMap_Area(Clustering_brainmask.infoCells.catAreaID(curIndCells(iCell)), :).*0.5, ...
            'MarkerFaceColor', cMap_Area(Clustering_brainmask.infoCells.catAreaID(curIndCells(iCell)), :), ...
            'LineWidth', 1) %cMap(iC,:), 
%         text(Y2(curIndCells(iCell),1)+1, Y2(curIndCells(iCell),2), catChanID(curIndCells(iCell),:))
        hold on
    end    
end
axis square
xlim([-50 40])
ylim([-25 25])
set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 2,  'Box', 'off')

% save
print(fig3a2, fullfile(dirFig, 'multipleFP_FigS_KMeansClustering_brainmaskVoxels_2DPCA_allCells_colorArea'), '-depsc')
print(fig3a2, fullfile(dirFig, 'multipleFP_FigS_KMeansClustering_brainmaskVoxels_2DPCA_allCells_colorArea'), '-r300', '-dtiff')

% Colored by cluster
fig3a3=figure;
set(fig3a3, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 700 700])
Y2(:,1) = score_e(:,1); %.*(-1); % to match it to the MDS plot
Y2(:,2) = score_e(:,2);
% Y3(:,1) = score_e(:,1); %.*(-1); % to match it to the MDS plot
% Y3(:,2) = score_e(:,2);
% Y3(:,3) = score_e(:,3);
for iC=1:curK
    
    iK = reorderCluster(iC);
    
    curIndCells = indSortChan_reorder(sortedClust_reorder==iK); %resultProbClustering(iK).validIndCells;
    combCells = nchoosek(curIndCells,2);

    figure(fig3a3);
%     line([Y2(combCells(:,1), 1) Y2(combCells(:,2),1)]', [Y2(combCells(:,1), 2) Y2(combCells(:,2), 2)]', 'Color', cMap(iC,:), 'LineWidth', 2);
    hold on;
    for iCell = 1:length(curIndCells)
        plot(Y2(curIndCells(iCell), 1), Y2(curIndCells(iCell), 2), 'o-', 'Marker', 'o', ... %marker{indMonkey(curIndCells(iCell))},...
            'MarkerSize', 10, 'MarkerEdgeColor', cMap(iC,:).*0.5, 'MarkerFaceColor', cMap(iC,:), 'LineWidth', 1) %
%         text(Y2(curIndCells(iCell),1)+1, Y2(curIndCells(iCell),2), catChanID(curIndCells(iCell),:))
        
%         plot3(Y3(curIndCells(iCell),1), Y3(curIndCells(iCell),2), Y3(curIndCells(iCell),3), 'o-',... %'Marker', marker{indMonkey(curIndCells(iCell))}, ...
%             'LineWidth', 2, 'MarkerSize', 12,...
%             'Color', cMap(iC,:), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(iC,:));
%         text(Y3(curIndCells(iCell),1)+1, Y3(curIndCells(iCell),2), Y3(curIndCells(iCell),3), catChanID(curIndCells(iCell),:))
        hold on
    end    
end
axis square
xlim([-50 40])
ylim([-25 25])
set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 2,  'Box', 'off')

% save
print(fig3a3, fullfile(dirFig, 'multipleFP_FigS_KMeansClustering_brainmaskVoxels_2DPCA_allCells_colorCluster'), '-depsc')
print(fig3a3, fullfile(dirFig, 'multipleFP_FigS_KMeansClustering_brainmaskVoxels_2DPCA_allCells_colorCluster'), '-r300', '-dtiff')



% %% Fig 3D: 2-D MDS plot showing K-means clustering results
% D = pdist(Clustering_brainmask.matR', 'euclidean');
% [Y2,stress,disparities] = mdscale(D,2);
% [Y3,stress,disparities] = mdscale(D,3);
% 
% curK=10; %7;
% locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
% locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
% [sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));
% 
% % matIndClust_SU = cat(2, Clustering_moviemask_valid.resultKMeans.SU_indCluster); %cat(2, Clustering.resultKMeans.SU_indCluster);
% % curK = 7; %
% % [sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
% indNewCluster = [4 2 3 6 7 5 1]; % [4 2 5 1 7 6 3]; %[1 2 3 4 5]; %[1 3 6 4 2 5]; %[2 4 1 3]; %[1 2 3 4 5]; %[3 4 2 5 1]; %[1 2 3 4 5 6 7]; %[4 1 6 3 5 2 7]; % cluster #4 is now cluster 1
% 
% % 
% 
% % catChanID = cat(1, paramClustering.validChanID);
% % indMonkey = str2num(catChanID(:,5));
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
% 
% % save
% print(fig3a2, fullfile(dirFig, 'figRev3D_2DMDS_diffMonk'), '-r300', '-dtiff')
% print(fig3a2, fullfile(dirFig, 'figRev3D_2DMDS_diffMonk'), '-depsc')
% 
% % fig3a2=figure;
% % set(fig3a2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 700 700])
% % for iC=1:curK
% %     
% %     iK = indNewCluster(iC);
% %     figure(fig3a2); 
% %     curChan = indSortChan(sortedClust==iK);
% %     plot3(Y3(curChan,1), Y3(curChan,2), Y3(curChan,3), 'o-','LineWidth', 2, 'MarkerSize', 10,...
% %         'MarkerEdgeColor','k', 'MarkerFaceColor', 'w', 'Color', cMap(iC,:)); %'MarkerFaceColor', cMap(iC,:), 'Color', cMap(iC,:));
% %     hold on;
% %     plot3(Y3(resultProbClustering(iC).validIndCells,1), Y3(resultProbClustering(iC).validIndCells,2), Y3(resultProbClustering(iC).validIndCells,3),...
% %         'o','LineWidth', 2, 'MarkerSize', 10,...
% %         'MarkerEdgeColor','k', 'MarkerFaceColor', cMap(iC,:), 'Color', cMap(iC,:));
% % %     text(Y2(curChan,1)+1, Y2(curChan,2), catChanID(curChan,:)) %paramCorr.validChanID(curChan,:))
% %     hold on;
% % end

%% Fig 3a_1: Correlation matrix between voxels and cells
% fig3a1=figure;
% set(fig3a1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 300 200])
% imagesc(matR_SU(1:11, 1:7)') %imagesc(matR_SU(1:20, 1:10)')
% axis off
% rgb=hslcolormap(256, 'bc.yr', 1, [0.2 1 0.2]); colormap(rgb)
% caxis([-0.18 .18])
% % save
% print(fig3a1, fullfile(dirFig, 'fig3a1'), '-r150', '-dtiff')
% 
% % % make blue-white-red colorbar
% % cval = 0.25;
% % cmin = -cval; cmax = cval;
% % colornum = 256;
% % colorInput = [1 0 0; 1 1 1; 0 0 1];
% % oldSteps = linspace(-1, 1, length(colorInput));
% % newSteps = linspace(-1, 1, colornum);
% % for j=1:3 % RGB
% %     newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1);
% % end
% % endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
% % newmap = squeeze(newmap_all(1:endPoint, :));
% % figure(gcf)
% % set(gca, 'CLim', [cmin cmax])
% % colormap(flipud(newmap))
% % set(gca, 'TickDir', 'out')
% % box off
% % c=colorbar;




%% 3D MDS
% [Y,stress,disparities] = mdscale(D,3);
%     
% matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster); %cat(2, Clustering.resultKMeans.SU_indCluster);
% curK = 5; %6; %7; %4;
% [sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
% 
% fig3a2=figure;
% set(fig3a2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1300 600 700 700])
% % Using PLOT3: main MDS space
% indNewCluster = [1 2 3 4 5 6 7]; %[4 1 6 3 5 2 7]; % cluster #4 is now cluster 1
% 
% for iC=1:curK
%     
%     iK = indNewCluster(iC);
%     figure(fig3a2); 
%     curChan = indSortChan(sortedClust==iK);
%     p_org=plot3(Y(curChan,1), Y(curChan,2), Y(curChan,3), 'o-','LineWidth', 2, 'MarkerSize', 12,...
%         'Color', cMap(iC,:), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(iC,:)); 
%     text(Y(curChan,1)+1, Y(curChan,2), Y(curChan,3), paramCorr.validChanID(curChan,:))
%     hold on;
% end
% 
% % Using PATCH
% % vert = Y;
% faces = padcat(indSortChan(sortedClust==1), indSortChan(sortedClust==2), indSortChan(sortedClust==3),...
%     indSortChan(sortedClust==4), indSortChan(sortedClust==5), indSortChan(sortedClust==6),...
%     indSortChan(sortedClust==7));
% faces = faces';
% % p_org = patch('Faces', faces, 'Vertices', vert, 'Marker', 'o'); hold on;
% % set(p_org, 'FaceColor', 'none', 'EdgeColor', 'flat', 'MarkerFaceColor', 'flat', 'FaceVertexCData', cdata)
% 
% % Draw projections of each dimension
% % set ideal view
% view([168.5 28])  %      view([44 44])
% % get axis limits
% axis tight
% xl = get(gca, 'XLim');
% yl = get(gca, 'YLim');
% zl = get(gca, 'ZLim');
% 
% % coordinates for background
% % x-y plane (bottom)
% bg_xy = [xl(2) xl(2) xl(1) xl(1); yl(1) yl(2) yl(2) yl(1); zl(1) zl(1) zl(1) zl(1)]; % background coords for xy plane
% vert_xy = cat(2, Y(:,1:2), repmat(zl(1), size(Y,1),1)); % data coords for xy plane
% % x-z plane (side)
% bg_xz = [xl(2) xl(2) xl(1) xl(1); yl(1) yl(1) yl(1) yl(1); zl(2) zl(1) zl(1) zl(2)]; % background coords for xy plane
% vert_xz = cat(2, Y(:,1), repmat(yl(2), size(Y,1),1), Y(:,3)); % data coords for xz plane
% % y-z plane (side)
% bg_yz = [xl(1) xl(1) xl(1) xl(1); yl(1) yl(1) yl(2) yl(2); zl(2) zl(1) zl(1) zl(2)]; % background coords for xy plane
% vert_yz = cat(2,repmat(xl(1), size(Y,1), 1), Y(:,2:3)); % data coords for yz plane
% 
% figure(fig3a2); hold on
% % grayVal = 0.95;
% % grayBG_xy = fill3(bg_xy(1,:), bg_xy(2,:), bg_xy(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); 
% % grayBG_xz = fill3(bg_xz(1,:), bg_xz(2,:), bg_xz(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); % for y-z plane
% % grayBG_yz = fill3(bg_yz(1,:), bg_yz(2,:), bg_yz(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); % for y-z plane
% p(1) = patch('Faces', faces, 'Vertices', vert_xy, 'Marker', 'o'); hold on;
% % p(2) = patch('Faces', faces, 'Vertices', vert_xz, 'Marker', 'o'); hold on;
% % p(3) = patch('Faces', faces, 'Vertices', vert_yz, 'Marker', 'o'); hold on;
% 
% indOldCluster = [2 6 4 1 5 3 7]; % kmeans clustering #1 is now cell group 2 in main figure
% cdata=cMap(indOldCluster(matIndClust_SU(:,curK-1)),:);
% set(p, 'FaceColor', 'none', 'EdgeColor', 'none', 'MarkerFaceColor', 'flat', 'FaceVertexCData', cdata, 'MarkerSize', 4)
% 
% % do some extra work to make it look better
% % line for y-z plane
% l=line([bg_yz(1,:); bg_yz(1,[2 3 4 1])], [bg_yz(2,:); bg_yz(2, [2 3 4 1])], [bg_yz(3,:); bg_yz(3, [2 3 4 1])], 'Color', 'k', 'LineWidth', 2);
% set(l(1:2), 'LineStyle', ':')
% % line for x-z plane
% l=line([bg_xz(1,:); bg_xz(1,[2 3 4 1])], [bg_xz(2,:); bg_xz(2, [2 3 4 1])], [bg_xz(3,:); bg_xz(3, [2 3 4 1])], 'Color', 'k', 'LineWidth', 2);
% set(l(2:3), 'LineStyle', ':')
% % line for x-y plane
% l=line([bg_xy(1,:); bg_xy(1,[2 3 4 1])], [bg_xy(2,:); bg_xy(2, [2 3 4 1])], [bg_xy(3,:); bg_xy(3, [2 3 4 1])], 'Color', 'k', 'LineWidth', 2);
% set(l(3:4), 'LineStyle', ':')
% 
% % xlabel('axis 1'); ylabel('axis 2'); zlabel('axis 3'); % just for check
% set(gca, 'XTick', [], 'YTick', [], 'ZTick', [], 'Box', 'off')


%% 2d MDS
% [Y2,stress,disparities] = mdscale(D,2);
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% plot(Y2(:,1),Y2(:,2),'o','LineWidth',2, 'MarkerSize', 8);
% text(Y2(:,1)+1, Y2(:,2), paramCorr.validChanID)
% % 
% figMDS=figure;
% set(figMDS, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 600 550 550])
% for iK = 1:8
%     
%     matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster); %
%     curK = paramClustering.setK(iK); %Clustering.setK(iK);
%     [sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
%     
%     for ii=1:curK
%         
%         figure(figMDS);
%         curChan = indSortChan(sortedClust==ii);
%         plot(Y2(curChan,1), Y2(curChan,2), 'o-','LineWidth', 2, 'MarkerSize', 8,...
%             'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:), 'Color', cMap(ii,:));
% %         text(Y2(curChan,1)+1, Y2(curChan,2), paramCorr.validChanID(curChan,:))
%         hold on;
%     end
%     
%     axis square
%     title(sprintf('MDS plot: cluster # = %d', curK))
%     input('')
%     
% end

% alpha = 0.7;
% pm = cat(1, p.MarkerHandle); % get hidden marker handle to control transparency of vertex markers
% faceCData = uint8(cMap(sortedClust,:).*255)';
% faceCData(4,:) = uint8(ones(1, size(faceCData,2)).*alpha*255);
% 
% [pm.FaceColorData] = deal(faceCData);
% [pm.Size] = deal(8);



% for ii=1:curK
%     
%     figure(figMDS3d);
%     curChan = indSortChan(sortedClust==ii);
%     p_org=plot3(Y(curChan,1), Y(curChan,2), Y(curChan,3), 'o-','LineWidth', 2, 'MarkerSize', 10,...
%         'Color', cMap(ii,:), 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(ii,:));
% %     text(Y(curChan,1)+1, Y(curChan,2), Y(curChan,3), paramCorr.validChanID(curChan,:))
%     hold on;
% end
% 
% % set a perfect view before you get the limit of each axis
% view([40 30])
% 
% % get limits
% xl = get(gca, 'XLim');
% yl = get(gca, 'YLim');
% zl = get(gca, 'ZLim');
% 
% yy = [yl(1) yl(1) yl(2) yl(2)];
% xx = [xl(2) xl(2) xl(2) xl(2)];
% zz = [zl(2) zl(1) zl(1) zl(2)];
% ff = fill3(xx, yy, zz, ones(1,3).*0.95, 'EdgeColor', 'none'); % for y-z plane
% hold on;


% %% using plot
% l = size(Y, 1);
%     
% % draw other projections
% alpha=0.5;
% set(gcf, 'Renderer', 'OpenGL')
% % x-y plane
% for ii=1:curK    
%     figure(figMDS3d);
%     curChan = indSortChan(sortedClust==ii);
%     p_xy=plot3(Y(curChan,1), Y(curChan,2), repmat(zl(1), length(curChan), 1), 'o','LineWidth', 2, 'MarkerSize', 8,...
%         'Color', cMap(ii,:), 'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:));
% %     p_xy.MarkerHandle.FaceColorData(4,:) = uint8(ones(1, size(p_xy.MarkerHandle.FaceColorData,2)).*alpha*255);
%     hold on;
% end
% % x-z plane
% for ii=1:curK    
%     figure(figMDS3d);
%     curChan = indSortChan(sortedClust==ii);
%     plot3(Y(curChan,1), repmat(yl(2), length(curChan), 1), Y(curChan,3), 'o-','LineWidth', 2, 'MarkerSize', 8,...
%         'Color', cMap(ii,:), 'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:));
%     hold on;
% end
% % y-z plane
% for ii=1:curK    
%     figure(figMDS3d);
%     curChan = indSortChan(sortedClust==ii);
%     p=plot3(repmat(xl(2), length(curChan), 1), Y(curChan,2), Y(curChan,3), 'o','LineWidth', 2, 'MarkerSize', 8,...
%         'Color', cMap(ii,:), 'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:));
%     hold on;
% end
%     
%     axis square
%     grid on
%     title(sprintf('MDS plot: cluster # = %d', curK))
% 





% figure(figMDS3d); hold on
% grayVal = 0.95;
% grayBG_xy = fill3(bg_xy(1,:), bg_xy(2,:), bg_xy(3,:), ones(1,3).*grayVal, 'EdgeColor', 'none'); 
% % x-y plane
% for ii=1:curK    
%     figure(figMDS3d);
%     curChan = indSortChan(sortedClust==ii);
%     p_xy=plot3(Y(curChan,1), Y(curChan,2), repmat(zl(1), length(curChan), 1), 'o','LineWidth', 2, 'MarkerSize', 8,...
%         'Color', cMap(ii,:), 'MarkerEdgeColor', cMap(ii,:), 'MarkerFaceColor', cMap(ii,:));
% %     p_xy.MarkerHandle.FaceColorData(4,:) = uint8(ones(1, size(p_xy.MarkerHandle.FaceColorData,2)).*alpha*255);
%     hold on;
% end




%     input('')
%     
% end

% for iK=1:7
% curChanInd = indSortChan(sortedClust==iK);
% plot3(Y(curChanInd,1),Y(curChanInd,2), Y(curChanInd,3),'o','LineWidth',2, 'MarkerSize', 8, ...
%     'MarkerEdgeColor', cMap(iK,:), 'MarkerFaceColor', cMap(iK,:)); %, ...
% %     'Marker', marker{iK});
% hold on
% % text(Y(curChanInd,1)+1, Y(curChanInd,2), Y(curChanInd,3), paramCorr.validChanID(curChanInd,:))
% end


%% Plotting parameters
% % cMap = [0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./255;
% cMap = [228	26	28;
% 55	126	184;
% 77	175	74;
% 152	78	163;
% 255	127	0;
% 255 217 47; % dark yellow %255	255	255; %white %255	255	51; % yellow was too similar to another yellow in mat2
% 166	86	40;
% 247	129	191;
% ]./255;
% marker = {'o', '^', 'v', '<', '>', 'square', 'diamond'}; %{'o', '^', 'square', 'diamond'}; % {'o','*', 'x', 's', 'd', '+', '^'};
% % orderIndNewClust = [3 2 7 6 5 1 4];
% % cMap_newclust = cMap(orderIndNewClust,:);
% 
% %% Load data
% % 1) Clustering results
% load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  
% % load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked.mat', nameSubjNeural, nameSubjBOLD)));
% % 2) Movie-driven mask
% load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');
% % 3) Valid cells considering likelihood being clustered together
% load(fullfile(dirDataNeural, 'Clustering_SU_allCells_validVoxels_critCorr1_7Means_prob.mat'))

