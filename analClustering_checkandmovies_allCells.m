% doClustering_checkandmovies_allCells.m

% 
addpath('/library/matlab_utils/')

% Load the data
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was'}; %{'Dav', 'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan'}; %{'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjNeural = 'Was'; % 'Dan'; %'Spi';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = '/procdata/parksh/_macaque';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Directory for saving figures as graphic files
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

%% Load data
% 1) Clustering based on corr maps
load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  
% load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  

% Clustering based on time series
load(fullfile(dirDataNeural, sprintf('ClusteringSDF_%sMovie123.mat', cell2mat(setNameSubjNeural))))
% 2) Movie-driven mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');
% 3) Valid cells considering likelihood being clustered together
load(fullfile(dirDataNeural, 'Clustering_SU_allCells_validVoxels_critCorr1_7Means_prob.mat'))


% % Clustering based on corr maps
% load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked.mat', nameSubjNeural, nameSubjBOLD)));
% % Clustering based on time series
% load(fullfile(dirDataNeural, sprintf('ClusteringSDF_%sMovie123.mat', nameSubjNeural)))
% 
% load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr')


%% 
setK = paramClustering_global.setK; %Clustering.setK;

matWSS_corrMap=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS_corrMap(:,iK) = sum(Clustering_moviemask_valid.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
    
%     matWSS_SDF(:,iK) = sum(ClusteringSDF.resultKMeans(iK).SU_sumD); %
%     matWSS_SDFnorm(:,iK) = sum(ClusteringSDFnorm.resultKMeans(iK).SU_sumD); %
%     matWSS_SDF_TR(:,iK) = sum(ClusteringSDF_TR.resultKMeans(iK).SU_sumD); %
%     matWSS_SDF_MION(:,iK) = sum(ClusteringSDF_MION.resultKMeans(iK).SU_sumD); %
end

totalSS_corrMap = Clustering_moviemask_valid.totalSS_SU;
% betweenSS_corrMap = totalSS_corrMap-matWSS_corrMap;
propExplained_corrMap = (totalSS_corrMap-matWSS_corrMap)./totalSS_corrMap; %matExpVar./totalSS;

% totalSS_SDF = ClusteringSDF.totalSS;
% propExplained_SDF = (totalSS_SDF-matWSS_SDF)./totalSS_SDF; %matExpVar./totalSS;
% 
% totalSS_SDFnorm = ClusteringSDFnorm.totalSS;
% propExplained_SDFnorm = (totalSS_SDFnorm-matWSS_SDFnorm)./totalSS_SDFnorm; %matExpVar./totalSS;
% 
% totalSS_SDF_TR = ClusteringSDF_TR.totalSS;
% propExplained_SDF_TR = (totalSS_SDF_TR-matWSS_SDF_TR)./totalSS_SDF_TR; %matExpVar./totalSS;

% totalSS_SDF_MION = ClusteringSDF_MION.totalSS;
% propExplained_SDF_MION = (totalSS_SDF_MION-matWSS_SDF_MION)./totalSS_SDF_MION; %matExpVar./totalSS;


%% Compare different clustering
K=4; %8; %7;
locMode_corrMap = find(propExplained_corrMap(:,K-1)==mode(propExplained_corrMap(:,K-1)));
% locMode_SDF = find(propExplained_SDF(:,K-1)==mode(propExplained_SDF(:,K-1)));
% locMode_SDFnorm = find(propExplained_SDFnorm(:,K-1)==mode(propExplained_SDFnorm(:,K-1)));
% locMode_SDF_TR = find(propExplained_SDF_TR(:,K-1)==mode(propExplained_SDF_TR(:,K-1)));
% locMode_SDF_MION = find(propExplained_SDF_MION(:,K-1)==mode(propExplained_SDF_MION(:,K-1)));

indClust_SU = Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode_corrMap(1)); % based on corr map
% indClustSDF_SU = ClusteringSDF.resultKMeans(K-1).SU_indCluster(:, locMode_SDF(1)); % SDF in 10 Hz
% indClustSDFnorm_SU = ClusteringSDFnorm.resultKMeans(K-1).SU_indCluster(:, locMode_SDFnorm(1)); % zscored SDF in 10 Hz
% indClustSDFTR_SU = ClusteringSDF_TR.resultKMeans(K-1).SU_indCluster(:, locMode_SDF_TR(1)); % SDF in TR resolution
% indClustSDFMION_SU = ClusteringSDF_MION.resultKMeans(K-1).SU_indCluster(:, locMode_SDF_MION(1)); % MION convolved SDF

sortTargetK = K; %5; % 7;
[sortedClust, indSortChan]=sort(indClust_SU);
% [sortedClustSDF, indSortChanSDF]=sort(indClustSDF_SU);
% [sortedClustSDFNorm, indSortChanSDFNorm]=sort(indClustSDFnorm_SU);
% [sortedClustSDFTR, indSortChanSDFTR]=sort(indClustSDFTR_SU);
% [sortedClustSDFMION, indSortChanSDFMION]=sort(indClustSDFMION_SU);

% [sortedClust, indSortChan] = sort(Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode_corrMap(1)));
% [sortedClustSDF, indSortChanSDF] = sort(ClusteringSDF.resultKMeans(K-1).SU_indCluster(:, locMode_SDF(1)));
% [sortedClustSDFNorm, indSortChanSDFNorm] = sort(ClusteringSDFnorm.resultKMeans(K-1).SU_indCluster(:, locMode_SDFnorm(1)));
% [sortedClustSDFTR, indSortChanSDFTR] = sort(ClusteringSDF_TR.resultKMeans(K-1).SU_indCluster(:, locMode_SDF_TR(1)));
% [sortedClustSDFMION, indSortChanSDFMION] = sort(ClusteringSDF_MION.resultKMeans(K-1).SU_indCluster(:, locMode_SDF_MION(1)));


oldIndCluster = [1 2 5 6 3 7 4]; % [3 4 5 2 1]; % [4 1 6 3 5 2 7]; 
indSortChan_new = [];
for iC = 1:K %7
    curC = oldIndCluster(iC);
    tempind = indSortChan(sortedClust==curC);
    indSortChan_new = cat(1, indSortChan_new, tempind);
end
sortedClust_new = indClust_SU(indSortChan_new);
% indSortChan_org = indSortChan;
% indSortChan = indSortChan_new;


catChanID = cat(1, paramClustering.validChanID);
% Comparison figure
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [680    61   601   914])
image(cat(2, indClust_SU(indSortChan_new), indClustSDF_SU(indSortChan_new), ...
    indClustSDFnorm_SU(indSortChan_new), indClustSDFTR_SU(indSortChan_new), ...
    indClustSDFMION_SU(indSortChan_new)))
colormap(lines)
% h=colorbar;
% set(h, 'YLim', [0 sortTargetK]+0.5)
% set(h, 'YTick', [])
set(gca, 'XTick', [])
set(gca, 'YTick', 1:135, 'YTickLabel', catChanID(indSortChan_new,:))
set(gca, 'FontSize', 12)
title('Results of K-means clustering of cells')
set(gca, 'XTick', 1:5, 'XTickLabel', {'Corr maps', 'Time Series (raw)', 'Time Series (zscore)', 'Time Series (TR)', 'Time Series (MION)'}); % {'Corr maps', 'Time Series (fine)', 'Time Series (zscore)', 'Time Series (TR)', 'Time Series (MION)'})
xlabel('Basis of clustering')
ylabel('Cells')
% 
% print(gcf, fullfile(dirFig, sprintf('kmeans_clusteredSU_sortedK%d_comparison_allCells', sortTargetK)), '-r150', '-dtiff');



%% Probability of the clustering 
K=7; %8; %6; %7;
locMode = find(propExplained_corrMap(:,K-1)==mode(propExplained_corrMap(:,K-1)));
locMin = find(propExplained_corrMap(:,K-1)==min(propExplained_corrMap(:,K-1)));
indClust = Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode(1)); %
[sortClust_validVoxel, sortCell_validVoxel] = sort(Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode(1)));

% orderROI = [4 2 5]; %[4 1 3 2 5]; %% 1: parafoveal:ventral stream / 2: face patches / 3: Eccentric visual cortex / 4: Central V1 / 5: MT
orderClust_Cell = 1:K; %[1 2 5 6 3 7 4]; %1:K; % [1 2 5 6 3 7 4]; %[4 2 5 1 7 6 3]; %

reOrdSortCell=[];
for iCC=1:K
    reOrdSortCell = cat(1, reOrdSortCell, find(indClust==orderClust_Cell(iCC)));
end

fig3c=figure;
set(fig3c, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 510 470])
imagesc(Clustering_moviemask_valid.resultKMeans(K-1).matProb(reOrdSortCell, reOrdSortCell)); % how probable
axis square
axis off
set(gca, 'CLim', [0 1])

% % save
% print(fig3c, fullfile(dirFig, sprintf('figRev3C_Likelihood_%dMeans', K)), '-depsc')
% print(fig3c, fullfile(dirFig, sprintf('figRev3C_Likelihood_%dMeans', K)), '-dtiff', '-r200')
% % print(fig3c, fullfile(dirFig, sprintf('figRevS4_Likelihood_%dMeans', K)), '-dtiff', '-r200')


%% Compute feature correlation with cells
% load the cocatenated time courses in different resolution
setNameSubjNeural = {'Dav', 'Spi', 'Mat', 'Dan'}; 
load(fullfile(dirDataNeural, sprintf('matSDF_%s_Movie123.mat', cell2mat(setNameSubjNeural))), 'matSDF')

% Concatenate SDFs 
% matFR_TR = cat(2, matSDF.matFR_TR);
% matNeuralRGR = cat(2, matSDF.matNeuralRGR);
matFR_SU = cat(2, matSDF.matFR_SU); % in 10hz
matFR4fps = resample(matFR_SU, 4, 10);
clear matSDF matFR_SU
% matFR_SU_norm = cat(2, matSDF.matFR_SU_norm);

% Compute the average correlation for each cluster
% K=7;
% locMode_corrMap = find(propExplained_corrMap(:,K-1)==mode(propExplained_corrMap(:,K-1)));
% indClust_SU = Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode_corrMap(1)); % based on corr map
% [sortedClust, indSortChan]=sort(indClust_SU);
% oldIndCluster = [1 2 5 6 3 7 4]; % [3 4 5 2 1]; % [4 1 6 3 5 2 7]; 
% indSortChan_new = [];
% for iC = 1:K %7
%     curC = oldIndCluster(iC);
%     tempind = indSortChan(sortedClust==curC);
%     indSortChan_new = cat(1, indSortChan_new, tempind);
% end
% sortedClust_new = indClust_SU(indSortChan_new);

meanFRCluster4fps = []; %steFRCluster4fps=[];
oldIndCluster = [1 2 5 6 3 7 4]; 
sortTargetK = 7;
for iK = 1:sortTargetK
    indClust = oldIndCluster(iK);
    tempMatFR=[];
    tempMatFR = matFR4fps(:,indSortChan(sortedClust==indClust));
    meanFRCluster4fps(:,iK) = mean(tempMatFR, 2);
%     steFRCluster4fps(:,iK) = std(matFR(:,tempMatFR, 0, 2)./length(indSortChan(sortedClust==iK));
end

flagSM = 1; % flag for compression and smoothing
setMovie = [1 2 3];
fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
% ttt=load('/procdata/parksh/MovieRegressors/dbtmMriReg.mat'); % Face scale regressor (in TR unit)
% scaleRGR = ttt.reg.xx(7,:)';

% full regressors
catRGRfull=[];matRGRfull=[];
for iMov=1:length(setMovie)
    m = setMovie(iMov);
    matCurRGR = fullRGR4fps(m).smoRegressors; %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
    catRGRfull = cat(1, catRGRfull, matCurRGR); % concatenation across movies
end
% scaleRGR_resampled = resample(scaleRGR, 2.4*100, 0.25*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGRfull = catRGRfull; %cat(2, catRGRfull, scaleRGR_resampled); %cat(2, matRGR, scaleRGR);
varnamesfull = fullRGR4fps(1).features; % cat(1, fullRGR4fps(1).features, {'Face size'});

% subset of regressors
indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
% matRGRvalid = matRGRfull(:,indValidRGR);
varnamesvalid = varnamesfull(indValidRGR);

% Compute correlation between neural TS and feature TS
R_ClusterMovieRGRfull=NaN(size(matRGRfull,2), size(meanFRCluster4fps,2));
R_SUmovieRGRfull=NaN(size(matRGRfull,2), size(matFR4fps,2));
for iRGR = 1:size(matRGRfull,2)
    r_c=[]; r_su=[];
    
%     % averaged TS in each cluster
%     r_c = corr(matRGRfull(:,iRGR), meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');
    % single unit TS
    r_su = corr(matRGRfull(:,iRGR), matFR4fps, 'rows', 'complete', 'type', 'Spearman');
    
%     R_ClusterMovieRGRfull(iRGR, :) = r_c;
    R_SUmovieRGRfull(iRGR, :) = r_su;    
    
end

R_SUmovieRGRvalid = R_SUmovieRGRfull(indValidRGR,:);
R_ClusterMovieRGRvalid = R_ClusterMovieRGRfull(indValidRGR,:);


%% Plot
% reorder the full regressors 
indReorderRGR = [1 2 9 10 6 7 8 3 11:18 4 5 21 27 28 20 29 30 22 31 32:36 23 24 25 26 19];
catChanID = cat(1, paramClustering.validChanID);

% 1. Correlation matrix between clusters and regressors
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 650 880])
imagesc(R_ClusterMovieRGRfull(indReorderRGR, :));
colorbar;
set(gca, 'CLim', [-.6 .6], 'FontSize', 12)
set(gca, 'YTick', 1:length(varnamesfull), 'YTickLabel', varnamesfull(indReorderRGR,:))
xlabel('Cell Group')
title('Correlation between cluster time series and feature time series')
% make blue-white-red colorbar
cval = 0.6;
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

print(gcf, fullfile(dirFig, 'Corr_ClusterMovieRGR_movie123_ClusteringCorrMapTorRhoSigSpiArt'), '-dtiff', '-r150')

% 2. Correlation matrix between SUs and regressors
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1600 600])
imagesc(R_SUmovieRGRfull(indReorderRGR, :)) % indSortChan_new));
xPos = find(abs(diff(sortedClust_new))>0)+0.5;
line([xPos xPos]', repmat([0; 37], 1, 6), 'Color', 'k') 
set(gca, 'CLim', [-.6 .6], 'FontSize', 12)
set(gca, 'XTick', 1:135, 'XTickLabel', [])
set(gca, 'YTick', 1:length(varnamesfull), 'YTickLabel', varnamesfull(indReorderRGR,:))
% make blue-white-red colorbar
cval = 0.6;
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
title('Correlation between SU time series and feature time series')

print(gcf, fullfile(dirFig, 'Corr_SUMovieRGR_movie123_ClusteringCorrMapTorRhoSigSpiArt'), '-dtiff', '-r150')








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
dirMovie = '/procdata/parksh/Stimulus/Movies/Rhesus';
videoObj1 = VideoReader(fullfile(dirMovie, 'Movie1.avi'));
videoObj2 = VideoReader(fullfile(dirMovie, 'Movie2.avi'));
videoObj3 = VideoReader(fullfile(dirMovie, 'Movie3.avi'));

% % Prepare neural responses
% % First, get the cell response in 30fps time resolution 
% FR_dT30 = createCellRegressor_indMov_discreteTime(dirDataNeural, cellstr(paramCorr.validChanID),...
%     setMovie, 1/30); % in 
% % concatenate across movies
% matFR=[];
% for iUnit = 1:size(FR_dT30,1)
%     tempFR = cat(1, FR_dT30(iUnit, :).mnFR);
%     matFR(:,iUnit) = tempFR;
% end
% matFR_zscore = zscore(matFR);
% 
% % Then 
% % 1. collect the scenes that drives average SDF within cluster across cells
% meanFRCluster=[];steFRCluster=[];
% for iK = 1:sortTargetK
% tempMatFR=[];
% tempMatFR = matFR(:,indSortChan(sortedClust==iK));
% meanFRCluster(:,iK) = mean(matFR(:,indSortChan(sortedClust==iK)), 2);
% steFRCluster(:,iK) = std(matFR(:,indSortChan(sortedClust==iK)), 0, 2)./length(indSortChan(sortedClust==iK));
% end
% meanFRCluster_zscore = zscore(meanFRCluster);


% %% Make videos
% for iClust = 1:7
%     
%     % Find when activity was high (i.e. exceeds certain criterion in z-score)
%     critHigh = 3; %1.5; %2; % in z-score
%     locHigh = find(meanFRCluster_zscore(:,iClust)>critHigh); 
% %     critHigh = -1; %3; %1.5; %2; % in z-score
% %     locHigh = find(meanFRCluster_zscore(:,iClust)<critHigh);
%     if length(locHigh) < 1
%         continue;
%     end
%     % [i,j] = ind2sub(size(meanFRCluster_zscore), locHigh);
%     % taxis = 1:900; % in second
%     
%     % Preallocate the movie structure for three-movie length
%     clear reverseCorrMov
%     reverseCorrMov(1:length(locHigh)) = struct('cdata',zeros(vidHeight,vidWidth, 3,'uint8'), 'colormap',[]);
%     
%     % Collect the scenes
%     countFrame = 0;
%     for iMovie = 1:3
%         validFrame_range = [(5*(iMovie-1)*60*30)+1, (5*iMovie*60*30)];
%         validLocHigh = [];
%         validLocHigh = locHigh(locHigh>validFrame_range(1) & locHigh<validFrame_range(2));
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
% %     % Size a figure based on the video's width and height.
% %     hf = figure;
% %     set(hf, 'position', [150 150 videoObj1.Width  videoObj1.Height])
% %     % Play back the movie once at the video's frame rate.
% %     movie(hf, reverseCorrMov, 1, 10); %videoObj1.FrameRate);
%     
%     % Write a movie
%     newIndCluster = [3 2 7 6 5 1 4];
%     fileName = sprintf('scenesMov123_Cluster%d_meanClustSDF_zCrit%d_OldCluster%d.avi', newIndCluster(iClust), critHigh, iClust);
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