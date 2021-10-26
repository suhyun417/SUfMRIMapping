% doClustering_checkandmovies_masked.m

% 
addpath('/library/matlab_utils/')

% Load the data
nameSubjNeural = 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Clustering based on corr maps
load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked.mat', nameSubjNeural, nameSubjBOLD)));
% Clustering based on time series
load(fullfile(dirDataNeural, sprintf('ClusteringSDF_%sMovie123.mat', nameSubjNeural)))

load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr')

% Directory for saving figures as graphic files
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';


%% Compare different clustering
% Sort out cells based on particular K-means clustering 
matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster); % based on corr map
matIndClustSDF_SU = cat(2, ClusteringSDF.resultKMeans.SU_indCluster); % SDF in 10 Hz
matIndClustSDFnorm_SU = cat(2, ClusteringSDFnorm.resultKMeans.SU_indCluster); % zscored SDF in 10 Hz
matIndClustSDFMION_SU = cat(2, ClusteringSDF_MION.resultKMeans.SU_indCluster); % MION convolved SDF
matIndClustSDFTR_SU = cat(2, ClusteringSDF_TR.resultKMeans.SU_indCluster); % SDF in TR resolution
matIndClustSDFRGR_SU = cat(2, ClusteringSDFRGR.resultKMeans.SU_indCluster);

sortTargetK = 5; % 7;
[sortedClust, indSortChan]=sort(matIndClust_SU(:,sortTargetK-1));
[sortedClustSDF, indSortChanSDF]=sort(matIndClustSDF_SU(:,sortTargetK-1));
[sortedClustSDFNorm, indSortChanSDFNorm]=sort(matIndClustSDFnorm_SU(:,sortTargetK-1));
[sortedClustSDFMION, indSortChanSDFMION]=sort(matIndClustSDFMION_SU(:,sortTargetK-1));
[sortedClustSDFTR, indSortChanSDFTR]=sort(matIndClustSDFTR_SU(:,sortTargetK-1));
[sortedClustSDFRGR, indSortChanSDFRGR]=sort(matIndClustSDFRGR_SU(:,sortTargetK-1)); %(:,5-1));


oldIndCluster = [3 4 5 2 1]; % [4 1 6 3 5 2 7]; 
indSortChan_new = [];
for iC = 1:sortTargetK %7
    curC = oldIndCluster(iC);
    tempind = indSortChan(sortedClust==curC);
    indSortChan_new = cat(1, indSortChan_new, tempind);
end
indSortChan_org = indSortChan;
indSortChan = indSortChan_new;

% Comparison figure
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [680    61   601   914])
image(cat(2, matIndClust_SU(indSortChan,sortTargetK-1),... % matIndClustSDF_SU(indSortChan, sortTargetK-1), ...
    matIndClustSDFnorm_SU(indSortChan, sortTargetK-1), matIndClustSDFTR_SU(indSortChan, sortTargetK-1), ...
    matIndClustSDFMION_SU(indSortChan, sortTargetK-1), matIndClustSDFRGR_SU(indSortChan, sortTargetK-1)))
colormap(lines)
% h=colorbar;
% set(h, 'YLim', [0 sortTargetK]+0.5)
% set(h, 'YTick', [])
set(gca, 'XTick', [])
set(gca, 'YTick', 1:48, 'YTickLabel', paramCorr.validChanID(indSortChan,:))
set(gca, 'FontSize', 12)
title('Results of K-means clustering of cells')
set(gca, 'XTick', 1:5, 'XTickLabel', {'Corr maps', 'Time Series (zscore)', 'Time Series (TR)', 'Time Series (MION)', 'Movie feature correlation'}); % {'Corr maps', 'Time Series (fine)', 'Time Series (zscore)', 'Time Series (TR)', 'Time Series (MION)'})
xlabel('Basis of clustering')
ylabel('Cells')
% 
print(gcf, fullfile(dirFig, sprintf('kmeans_clusteredSU_sortedK%d_comparison_new_masked', sortTargetK)), '-r150', '-dtiff');


setKplot = sortTargetK-2:sortTargetK+2;
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [680    61   601   914])
image(matIndClust_SU(indSortChan,setKplot-1)) %matIndClustSDF_SU(indSortChanSDF,setKplot-1))
colormap(lines)
h=colorbar;
set(h, 'YLim', [0 setKplot(end)]+0.5)
set(h, 'YTick', [])
set(gca, 'YTick', 1:48, 'YTickLabel', paramCorr.validChanID(indSortChan,:)) %paramCorr.validChanID(indSortChanSDFRGR,:))
set(gca, 'FontSize', 12)
title('Results of K-means clustering of cells')
set(gca, 'XTick', 1:length(setKplot), 'XTickLabel', setKplot); %setKplot)
xlabel('Number of clusters')
ylabel('Cells')
print(gcf, fullfile(dirFig, sprintf('kmeans_clusteredSU_sortedK%d_newCorrMaps_masked', sortTargetK)), '-r150', '-dtiff');


% check the assignment across different k
matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster); % based on corr map
for sortTargetK = 4:6
[sortedClust, indSortChan]=sort(matIndClust_SU(:,sortTargetK-1));

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [680    61   100   914])
image(matIndClust_SU(indSortChan,sortTargetK-1)) %matIndClustSDF_SU(indSortChanSDF,setKplot-1))
colormap(lines)
set(gca, 'YTick', 1:48, 'YTickLabel', paramCorr.validChanID(indSortChan,:)) %paramCorr.validChanID(indSortChanSDFRGR,:))
set(gca, 'FontSize', 12)
end
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% image(matIndClustSDF_SU(indSortChanSDF,[2:7]-1))
% colormap(lines)
% h=colorbar;
% set(h, 'YLim', [0 setKplot(end)]+0.5)
% set(h, 'YTick', [])
% set(gca, 'XTick', [])
% set(gca, 'YTick', 1:48, 'YTickLabel', paramCorr.validChanID(indSortChanSDF,:))
% set(gca, 'FontSize', 12)
% title('Results of K-means clustering of cells based on timeseries')
% set(gca, 'XTick', 1:length([2:7]), 'XTickLabel', 2:7); %setKplot)
% xlabel('Number of clusters')
% ylabel('Cells')
% % print(gcf, fullfile(dirFig, sprintf('kmeans_clusteredSUSDF_sortedK%d', sortTargetK)), '-r150', '-dtiff');


%% Single unit time series
setMovie = [1 2 3];
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr') % we only need ID of valid channels

% First, get the cell response in 4fps time resolution to compare it with movie regressors
FR_dT = createCellRegressor_indMov_discreteTime(dirDataNeural, cellstr(paramCorr.validChanID),...
    setMovie, 0.25); % in 
% concatenate across movies
matFR=[];
for iUnit = 1:size(FR_dT,1)
    tempFR = cat(1, FR_dT(iUnit, :).mnFR);
    matFR(:,iUnit) = tempFR;
end

%
figure; % normalized time series sorted based on clustering result
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(zscore(matFR(:,indSortChan))')
set(gca, 'YTick', 1:48, 'YTickLabel', paramCorr.validChanID(indSortChan,:))
set(gca, 'XTick', 400:400:3600, 'XTickLabel', 100:100:900)
set(gca, 'FontSize', 12)
colormap(hot)
title('Normalized neural responses to movie 1 2 3')
ylabel('Cells')
xlabel('Time (s)')


%% Movie regressors
setMovie = [1 2 3];
flagSM = 1; % flag for compression and smoothing

fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
ttt=load('/procdata/parksh/MovieRegressors/dbtmMriReg.mat'); % Face scale regressor (in TR unit)
scaleRGR = ttt.reg.xx(7,:)';


% full regressors
catRGRfull=[];matRGRfull=[];
for iMov=1:length(setMovie)
    m = setMovie(iMov);
    matCurRGR = fullRGR4fps(m).smoRegressors; %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
    catRGRfull = cat(1, catRGRfull, matCurRGR); % concatenation across movies
end
scaleRGR_resampled = resample(scaleRGR, 2.4*100, 0.25*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGRfull = cat(2, catRGRfull, scaleRGR_resampled); %cat(2, matRGR, scaleRGR);
varnamesfull = cat(1, fullRGR4fps(1).features, {'Face size'});

% subset of regressors
indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
matRGRvalid = matRGRfull(:,indValidRGR);
varnamesvalid = varnamesfull(indValidRGR);



% reorder the full regressors 
indReorderRGR = [1 2 9 10 6 7 8 3 11:18 4 5 21 27 28 20 29 30 32 22 31 23 24 25 26 19];
% 'Luminance' 
%     'Contrast'
%     'Beta Contrast'
%     'Gamma Contrast'
%     'SF (0.0 to 0.2)'
%     'SF (3.0 to 100.0)'
%     'SF Ratio'
%     'Speed'
%     'Motion Contrast'
%     'Motion Dir Beta'
%     'Motion Dir Gamma'
%     'Motion Div Mean'
%     'Motion Div Rectify'
%     'Motion Div STD'
%     'Motion Div Beta'
%     'Motion Div Gamma'
%     'Rightward Motion'
%     'Leftward Motion'
%     'One face'
%     'One Face (full)'
%     'One Face (side view)'
%     'Faces'
%     'Faces (full)'
%     'Faces (side view)'
%     'Face size'
%     'Body parts'
%     'Hands'
%     'Conspecifics'
%     'Humans'
%     'Any animal'
%     'Dyadic interaction'
%     'Scene Cuts'
matRGRfull = matRGRfull(:, indReorderRGR); %cat(2, matRGR, scaleRGR);
varnamesfull = varnamesfull(indReorderRGR);


%% Cells and movie features
% Correlation between each cell and movie RGRs
[R_SUmovieRGR] = corr(matRGRfull, matFR, 'rows', 'complete');
[R_SUmovieRGRvalid] = corr(matRGRvalid, matFR, 'rows', 'complete');

% correlation matrix between SUs and regressors
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 800 800])
imagesc(R_SUmovieRGR');
set(gca, 'XTick', 1:length(fullRGR4fps(1).features)+1, 'XTickLabel', []) %varnamesfull)
set(gca, 'YTick', 1:length(paramCorr.validChanIndex), 'YTickLabel', cellstr(paramCorr.validChanID))

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 400 800])
imagesc(R_SUmovieRGRvalid(:, indSortChan)');
set(gca, 'XTick', 1:size(R_SUmovieRGRvalid,1), 'XTickLabel', []) %varnamesfull)
set(gca, 'YTick', 1:length(paramCorr.validChanIndex), 'YTickLabel', cellstr(paramCorr.validChanID(indSortChan,:)))

% make blue-white-red colorbar
cval = 0.7;
cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end
endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));
figure(gcf)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;
set(gca, 'FontSize', 12)
print(gcf, fullfile(dirFig, 'Corr_SUMovieRGR_movie123Tor_rgrReordered'), '-dtiff', '-r150')
print(gcf, fullfile(dirFig, 'Corr_SUMovieRGR_movie123Tor_rgrSubset'), '-dtiff', '-r150')

% correlation matrix between SUs and regressors
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 800 800])
imagesc(R_SUmovieRGR(:, indSortChan)');
set(gca, 'XTick', 1:length(fullRGR4fps(1).features)+1, 'XTickLabel', []) %varnamesfull)
set(gca, 'YTick', 1:length(paramCorr.validChanIndex), 'YTickLabel', cellstr(paramCorr.validChanID(indSortChan,:)))

% make blue-white-red colorbar
cval = 0.7;
cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end
endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));
figure(gcf)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;
set(gca, 'FontSize', 12)
print(gcf, fullfile(dirFig, 'Corr_SUMovieRGR_movie123Tor_ClusteringCorrMapTorArt'), '-dtiff', '-r150')

%% Correlation between cluster & movie regressor
% average across cells for each cluster
% matFR_zscored = zscore(matFR);
meanFRCluster=[];steFRCluster=[];
for iK = 1:sortTargetK
    tempMatFR=[];
    tempMatFR = matFR(:,indSortChan(sortedClust==iK));
    meanFRCluster(:,iK) = mean(matFR(:,indSortChan(sortedClust==iK)), 2);
    steFRCluster(:,iK) = std(matFR(:,indSortChan(sortedClust==iK)), 0, 2)./length(indSortChan(sortedClust==iK));
end

[R_avgClusterSU_RGR] = corr(matRGRfull, meanFRCluster, 'rows', 'complete'); %corr(matRGR, meanFRCluster, 'rows', 'complete');

% varnames_fullRGR = fullRGR4fps(1).features; %cat(1, fullRGR4fps(1).features, {'Face size'}); %cat(1, fullRGR4fps(1).features(indValidRGR), {'Face size'});

figCorrClusterMovie = figure;
set(figCorrClusterMovie, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(R_avgClusterSU_RGR);
set(gca, 'YTick', 1:size(matRGRfull,2), 'YTickLabel', varnamesfull);
set(gca, 'FontSize', 15)
% set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', fullRGR4fps(1).features(indValidRGR));
xlabel('Clusters')
ylabel('Features')
title('Correlation between averaged cluster timeseries and movie rgrs')

% make blue-white-red colorbar
cval = 0.8;
cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end
endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));

figure(figCorrClusterMovie)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;

% print(gcf, fullfile(dirFig, 'Corr_ClusterMovieRGR_movie123Tor'), '-dtiff', '-r150')



%% Clustering based on correlation between movie regressors
% setK = 2:15;
% opts = statset('Display','final');
% numReplicates = 5;
% 
% % matCluster_SU = NaN(size(matR_SU, 2), 
% for iK = 1:length(setK)
%     
%     K = setK(iK);
%     % Cluster single units based on correlation with movie regressors
%     [IDX_SUrgr, C, SUMD_SUrgr] = kmeans(R_SUmovieRGR', K, 'Replicates', numReplicates, 'Options', opts); 
%     % Cluster regressors based on 50 singel unit correlation
%     [IDX_rgr, C, SUMD_rgr] = kmeans(R_SUmovieRGR, K, 'Replicates', numReplicates, 'Options', opts);
%     
%     
%     ClusteringSDFRGR.resultKMeans(iK).SU_indCluster = IDX_SUrgr;
%     ClusteringSDFRGR.resultKMeans(iK).SU_sumD = SUMD_SUrgr;
%     ClusteringSDFRGR.resultKMeans(iK).rgr_indCluster = IDX_rgr;
%     ClusteringSDFRGR.resultKMeans(iK).rgr_sumD = SUMD_rgr;
%     
% end
% plot Clustering results
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% for iK=1:length(setK)
%     plot(iK+1, ClusteringSDFRGR.resultKMeans(iK).SU_sumD, 'ko', 'MarkerSize', 10, 'LineWidth', 2)
%     hold on
% end
% xlim([1 16])
% set(gca, 'LineWidth', 2, 'FontSize', 12)
% box off
% title('Clustering of single units based on movie regressor correlation')
% xlabel('Number of cluster')
% ylabel('Within-cluster distance')
% print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_SUSDFRGR'), '-depsc');
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% for iK=1:length(setK)
%     plot(iK+1, ClusteringSDFRGR.resultKMeans(iK).rgr_sumD, 'ko', 'MarkerSize', 10, 'LineWidth', 2)
%     hold on
% end
% xlim([1 16])
% set(gca, 'LineWidth', 2, 'FontSize', 12)
% box off
% title('Clustering of regressors')
% xlabel('Number of cluster')
% ylabel('Within-cluster distance')
% print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_rgr'), '-depsc');


%% TEMP: regressor maps
% catRGRTR=[];matRGRTR=[];
% for iMov=1:length(setMovID)
%     m = setMovID(iMov);
%     matCurRGRTR = resample(fullRGR4fps(1).smoRegressors,  0.25*100, 2.4*100); %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
%     catRGRTR = cat(1, catRGRTR, matCurRGRTR); % concatenation across movies
% end
% matRGRTR = catRGRTR; % cat(2, catRGRTR, scaleRGR); %cat(2, matRGR, scaleRGR);
% 
% 
% matR_RGR = NaN(nVox, size(matRGRTR,2));
% matRGRTR = matRGRTR-repmat(nanmean(matRGRTR), nt, 1); % centering
% matRGRTR = doConv(matRGRTR,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
% 
% 
% [Rvals, Pvals] = corr(reshape(fmritc, nVox, nt)', matRGRTR',...
%     'rows','complete', 'type', 'Spearman');
% 
% matR_RGR = Rvals.*(-1); % because of MION
% % mapR = reshape(Rvals, [nx, ny, nz]).*(-1); % because of MION
% %     

% setK = 2:15;
% opts = statset('Display','final');
% numReplicates = 5;
% 
% % matCluster_SU = NaN(size(matR_SU, 2), 
% for iK = 1:length(setK)
%     
%     K = setK(iK);
%     % Cluster regressors based on whole brain corr map
%     [IDX_rgr, C, SUMD_rgr] = kmeans(matR_RGR', K, 'Replicates', numReplicates, 'Options', opts);
%     
%     
%     ClusteringRGR_map.resultKMeans(iK).rgr_indCluster = IDX_rgr;
%     ClusteringRGR_map.resultKMeans(iK).rgr_sumD = SUMD_rgr;
% %     ClusteringSDFRGR.resultKMeans(iK).rgr_indCluster = IDX_rgr;
% %     ClusteringSDFRGR.resultKMeans(iK).rgr_sumD = SUMD_rgr;
%     
% end
% 
% % Sort out cells based on particular K-means clustering 
% matIndClustRGRmap= cat(2, ClusteringRGR_map.resultKMeans.rgr_indCluster);
% 
% sortTargetK = 7;
% [sortedClustRGRmap, indSortRGRmap]=sort(matIndClustRGRmap(:,sortTargetK-1));
% 
% setKplot = sortTargetK-2:sortTargetK+2;
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% image(matIndClustRGRmap(indSortRGRmap,setKplot-1))
% colormap(lines)
% h=colorbar;
% set(h, 'YLim', [0 setKplot(end)]+0.5)
% set(h, 'YTick', [])
% set(gca, 'YTick', 1:31, 'YTickLabel', varnamesfull(indSortRGRmap))
% set(gca, 'FontSize', 12)
% title('Results of K-means clustering of cells based on movie rgrs')
% set(gca, 'XTick', 1:length(setKplot), 'XTickLabel', setKplot); %setKplot)
% xlabel('Number of clusters')
% ylabel('Cells')
% print(gcf, fullfile(dirFig, sprintf('kmeans_clusteredSUSDF_sortedK%d', sortTargetK)), '-r150', '-dtiff');



%% Cluster and movie scenes
% Read videos (takes long time)
dirMovie = '/procdata/parksh/Movies';
videoObj1 = VideoReader(fullfile(dirMovie, 'Movie1.avi'));
videoObj2 = VideoReader(fullfile(dirMovie, 'Movie2.avi'));
videoObj3 = VideoReader(fullfile(dirMovie, 'Movie3.avi'));

% Prepare neural responses
% First, get the cell response in 30fps time resolution 
FR_dT30 = createCellRegressor_indMov_discreteTime(dirDataNeural, cellstr(paramCorr.validChanID),...
    setMovie, 1/30); % in 
% concatenate across movies
matFR=[];
for iUnit = 1:size(FR_dT30,1)
    tempFR = cat(1, FR_dT30(iUnit, :).mnFR);
    matFR(:,iUnit) = tempFR;
end
matFR_zscore = zscore(matFR);

% Then 
% 1. collect the scenes that drives average SDF within cluster across cells
meanFRCluster=[];steFRCluster=[];
for iK = 1:sortTargetK
tempMatFR=[];
tempMatFR = matFR(:,indSortChan(sortedClust==iK));
meanFRCluster(:,iK) = mean(matFR(:,indSortChan(sortedClust==iK)), 2);
steFRCluster(:,iK) = std(matFR(:,indSortChan(sortedClust==iK)), 0, 2)./length(indSortChan(sortedClust==iK));
end
meanFRCluster_zscore = zscore(meanFRCluster);


%% Make videos
for iClust = 1:7
    
    % Find when activity was high (i.e. exceeds certain criterion in z-score)
    critHigh = 3; %1.5; %2; % in z-score
    locHigh = find(meanFRCluster_zscore(:,iClust)>critHigh); 
%     critHigh = -1; %3; %1.5; %2; % in z-score
%     locHigh = find(meanFRCluster_zscore(:,iClust)<critHigh);
    if length(locHigh) < 1
        continue;
    end
    % [i,j] = ind2sub(size(meanFRCluster_zscore), locHigh);
    % taxis = 1:900; % in second
    
    % Preallocate the movie structure for three-movie length
    clear reverseCorrMov
    reverseCorrMov(1:length(locHigh)) = struct('cdata',zeros(vidHeight,vidWidth, 3,'uint8'), 'colormap',[]);
    
    % Collect the scenes
    countFrame = 0;
    for iMovie = 1:3
        validFrame_range = [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
        validLocHigh = [];
        validLocHigh = locHigh(locHigh>validFrame_range(1) & locHigh<validFrame_range(2));
        
        switch iMovie
            case 1
                obj = videoObj1;
            case 2
                obj = videoObj2;
            case 3
                obj = videoObj3;
        end
        
        % Read one frame at a time.
        for k = 1 : length(validLocHigh)
            indFrame = validLocHigh(k)-(5*(iMovie-1)*60*30);
            reverseCorrMov(countFrame+k).cdata = read(obj,indFrame);
        end
        countFrame = countFrame + length(validLocHigh);
    end
    
%     % Size a figure based on the video's width and height.
%     hf = figure;
%     set(hf, 'position', [150 150 videoObj1.Width  videoObj1.Height])
%     % Play back the movie once at the video's frame rate.
%     movie(hf, reverseCorrMov, 1, 10); %videoObj1.FrameRate);
    
    % Write a movie
    newIndCluster = [3 2 7 6 5 1 4];
    fileName = sprintf('scenesMov123_Cluster%d_meanClustSDF_zCrit%d_OldCluster%d.avi', newIndCluster(iClust), critHigh, iClust);
%     fileName = sprintf('scenesMov123_Cluster%d_meanClustSDF_zCritLow%d_OldCluster%d.avi', newIndCluster(iClust), critHigh, iClust);
    myObj = VideoWriter(fileName);
    myObj.FrameRate = 10;
%     myObj.Path = dirMovie;
    open(myObj);
    writeVideo(myObj, reverseCorrMov)
    close(myObj);
end

% Move files to /procdata
movefile('./scenesMov123*.avi', dirMovie)

%% Data mining: compare time series & movie scenes, instead of selecting out scenes
% Get the scene info from DM's results
for iMovie = 2:3
    
    switch iMovie
        case 1
            obj = videoObj1;
        case 2
            obj = videoObj2;
        case 3
            obj = videoObj3;
    end
    
    load(sprintf('/procdata/parksh/MovieRegressors/annotationMovie%d.mat', iMovie))
    sceneInfo = cat(2, sta', sto');
    nScene = size(sceneInfo, 1);
    
    oldOrderCluster = [6 2 1 7 5 4 3];
    
    figTS = figure;
    set(figTS, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 900])
    SP(1) = subplot('Position', [0.1300    0.65    0.7750    0.3]); %[0.1 0.75 0.8 0.2]); % for time series
    SP(2) = subplot('Position', [0.1 0.05 0.5 0.5]); % for movie image
    
    % % Size a figure based on the video's width and height.
    % hf = figure;
    % set(hf, 'Color', 'w','position', [150 150 obj.Width  obj.Height])
    
    for iScene = 1:nScene
        
        setFrames = sceneInfo(iScene,1):sceneInfo(iScene,2);
        
        
        subplot(SP(1))
        plot(SP(1), setFrames, meanFRCluster_zscore(setFrames+(5*(iMovie-1)*60*30), oldOrderCluster), 'o-')
        LG= legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5', 'Cluster 6', 'Cluster 7',...
            'Location', [0.0056    0.8433    0.0804    0.1497]);
        axis tight
        ylim([-1 10])
        xlabel('Frame #')
        ylabel('Normalized response (z)')
        title(sprintf('Movie %d: Scene %d', iMovie, iScene))
        
        frame = read(obj,setFrames(1));
        a(1).cdata = frame;
        a(1).colormap = [];
        imageFrame = frame2im(a);
        subplot(SP(2));
        imagesc(imageFrame);
        
        while 1
            [x, y, but] = ginput(1);
            indFrame = round(x);
            frame = read(obj,indFrame);
            a(1).cdata = frame;
            a(1).colormap = [];
            imageFrame = frame2im(a);
            
            %         figure(hf);
            subplot(SP(2));
            imagesc(imageFrame);
            
            if but ~= 1
                break;
            end
            
        end
    end
end


% 
% % 2. take the scenes drive all the cells in the cluster
% for iClust = 1:7
%     
%     matFR_zscore_cluster=[];
%     matFR_zscore_cluster = matFR_zscore(:,indSortChan(sortedClust==iClust));
%     
%     % Find when activity exceeds certain criterion in z-score in ALL the
%     % cells in this cluster (i.e. takes intersection from all of the cells)
%     critHigh = 3; %1.5; %2; % in z-score
%     locHigh = matFR_zscore_cluster>critHigh; 
%     intersectScene = find(sum(locHigh,2)>2);
% % %     critHigh = -1; %3; %1.5; %2; % in z-score
% % %     locHigh = find(meanFRCluster_zscore(:,iClust)<critHigh);
%     if length(intersectScene) < 1
%         continue;
%     end
% %     % [i,j] = ind2sub(size(meanFRCluster_zscore), locHigh);
% %     % taxis = 1:900; % in second
%     
%     % Preallocate the movie structure for three-movie length
%     clear reverseCorrMov
%     reverseCorrMov(1:length(intersectScene)) = struct('cdata',zeros(vidHeight,vidWidth, 3,'uint8'), 'colormap',[]);
%     
%     % Collect the scenes
%     countFrame = 0;
%     for iMovie = 1:3
%         validFrame_range = [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
%         validLocHigh = [];
%         validLocHigh = intersectScene(intersectScene>validFrame_range(1) & intersectScene<validFrame_range(2));
%         
%         switch iMovie
%             case 1
%                 obj = videoObj1;
%             case 2
%                 obj = videoObj2;
%             case 3
%                 obj = videoObj3;
%         end
%         
%         % Read one frame at a time.
%         for k = 1 : length(validLocHigh)
%             indFrame = validLocHigh(k)-(5*(iMovie-1)*60*30);
%             reverseCorrMov(countFrame+k).cdata = read(obj,indFrame);
%         end
%         countFrame = countFrame + length(validLocHigh);
%     end
%     
%     % Size a figure based on the video's width and height.
%     hf = figure;
%     set(hf, 'position', [150 150 videoObj1.Width  videoObj1.Height])
%     % Play back the movie once at the video's frame rate.
%     movie(hf, reverseCorrMov, 1, 10); %videoObj1.FrameRate);
%     
%     % Write a movie
%     newIndCluster = [3 2 7 6 5 1 4];
%     fileName = sprintf('scenesMov123_intersectClust_Cluster%d_zCrit%d_OldCluster%d.avi', newIndCluster(iClust), critHigh, iClust);
% %     fileName = sprintf('scenesMov123_Cluster%d_meanClustSDF_zCritLow%d_OldCluster%d.avi', newIndCluster(iClust), critHigh, iClust);
%     myObj = VideoWriter(fileName);
%     myObj.FrameRate = 10;
% %     myObj.Path = dirMovie;
%     open(myObj);
%     writeVideo(myObj, reverseCorrMov)
%     close(myObj);
% end
% 
% % Move files to /procdata
% movefile('./scenesMov123*.avi', dirMovie)


% % Visualize the clustered voxels
% indClustVox_ConvR = -0.6+(matIndClust_Vox(:,4)-1).*(0.4);
% mapClustVox = reshape(indClustVox_ConvR, [nx, ny, nz]);

%% Some visualization fun
% Plot clustering result on arbitrarily assigned axis
% first, set the coordinates for each cell: main cluster is whole
% brain corr clustering
% center for each cluster
centerCoords  = [0 0;...
    0 1;...
    1 0;...
    0 -1;...
    -1 0;...
    0 2;...
    2 0];
centerCoords = centerCoords.*1.5;

% theta = linspace(0,2*pi,150);
% x = sin(theta) + 0.75*rand(1,150);
% y = cos(theta) + 0.75*rand(1,150);

xCoord = []; yCoord = [];
for iK = 1:7
    tempx=[]; tempy=[];
    
%     theta = linspace(0,2*pi,length(find(sortedClust==iK)))';
%     tempx = centerCoords(iK,1) + sin(theta) + 0.75.*rand(length(find(sortedClust==iK)),1);
%     tempy = centerCoords(iK,2) + cos(theta) + 0.75.*rand(length(find(sortedClust==iK)),1);

    tempx = centerCoords(iK,1) + 0.5.*randn(length(find(sortedClust==iK)),1);
    tempy = centerCoords(iK,2) + 0.5.*randn(length(find(sortedClust==iK)),1);
    xCoord = cat(1, xCoord, tempx);
    yCoord = cat(1, yCoord, tempy);
end

% Clustering based on whole brain
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [676   250   586   583])
scatter(xCoord, yCoord, 100, sortedClust, 'LineWidth', 3)
colormap(jet(15))
axis square
set(gca, 'XTick', [], 'YTick', [])
set(gca, 'FontSize', 15)
title('Clustering based on whole brain correlation map')
%  print(gcf, fullfile(dirFig, 'kmeans_scatterplot_wholebrain'), '-depsc');

 % CLustering based on SDF
 [lia1, locb1] = ismember(indSortChan, indSortChanSDF);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [676   250   586   583])
scatter(xCoord, yCoord, 100, sortedClustSDF(locb1), 'LineWidth', 3)
colormap(jet(15))
axis square
set(gca, 'XTick', [], 'YTick', [])
set(gca, 'FontSize', 15)
title('Clustering based on cell time series')
%  print(gcf, fullfile(dirFig, 'kmeans_scatterplot_SDF'), '-depsc');
 
 
  % CLustering based on movie RGRs
  [lia1, locb2] = ismember(indSortChan, indSortChanSDFRGR);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [676   250   586   583])
scatter(xCoord, yCoord, 100, sortedClustSDFRGR(locb2), 'LineWidth', 3)
colormap(jet(15))
axis square
set(gca, 'XTick', [], 'YTick', [])
set(gca, 'FontSize', 15)
title('Clustering based on correlation with movie rgrs')
%  print(gcf, fullfile(dirFig, 'kmeans_scatterplot_movieRGR'), '-depsc');