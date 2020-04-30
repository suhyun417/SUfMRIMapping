% doClusteringSDF.m
%

% load the data
nameSubjNeural = 'Spi'; %'Tor';
% nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
% 
addpath('/library/matlab_utils/')

% Load data files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
% dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% fileNameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
% filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
load(fullfile(dirDataNeural, filenameNeural))
% load(fullfile(dirDataNeural, fileNameNeural_BLP))
% fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
% load(fullfile(dirDataBOLD, filenameBOLD))

% % Get movie IDs common in two dataset
% commonSetMovie = intersect(paramBOLD.unimov, paramSDF.setMovIDs);
% dataBOLD.mvoltc = voltcIndMov(commonSetMovie);
% dataBOLD.unimov = commonSetMovie;

setMovie = [1 2 3];
switch lower(nameSubjNeural)
    case 'spi'
        excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
        validC = setdiff(1:length(paramSDF.setCellIDs), excChanIndex)';        
    otherwise
        [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
        validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
end 
indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);


% MION function
TR=2.4;
k = gampdf([-40:TR:40],4,2);

% % 1. fMRI tc
% fmritc=[];
% indMovieBOLD = find(ismember(dataBOLD.unimov, setMovie)>0);
% for iM = indMovieBOLD %1:length(indMovieBOLD)
%     curvoltc = dataBOLD.mvoltc{iM};
%     avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
%     if ~isempty(find(avgvoltc==0, 1))
%         avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
%     end
%     pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
%     fmritc = cat(4,fmritc,pcvoltc);
% end

%% 1-1. K-means Clustering based on timeseries 
% (1) MION-convolved timeseries (neuralrgrs), i.e. in TR resolution
nt = 375;
matNeuralRGR = NaN(nt, length(validC));
for iChan = 1:length(validC) % compute correlation channel-by-channel
    neuralrgrs=[];
    neuralrgrs = cat(1, S(validC(iChan), indMovieNeuron).mnFR);
    neuralrgrs = neuralrgrs-mean(neuralrgrs); % centering
    neuralrgrs = doConv(neuralrgrs,k); % convolve MI
    matNeuralRGR(:,iChan) = neuralrgrs';
end

% no MION
matFR_TR = NaN(nt, length(validC));
for iChan = 1:length(validC) % compute correlation channel-by-channel
    tempFR = [];
    tempFR = cat(1, S(validC(iChan), indMovieNeuron).mnFR);    
    matFR_TR(:,iChan) = tempFR;
end

setK = 2:15;
opts = statset('Display','final');
numReplicates = 100; %5;


% matCluster_SU = NaN(size(matR_SU, 2), 
for iK = 1:length(setK)
    
    K = setK(iK);
    % Cluster single units based on time series
    [IDX_SU, C, SUMD_SU] = kmeans(matNeuralRGR', K, 'Replicates', numReplicates, 'Options', opts); 
    % Cluster time points based on 50 singel unit correlation
    [IDX_t, C, SUMD_t] = kmeans(matNeuralRGR, K, 'Replicates', numReplicates, 'Options', opts);
        
    ClusteringSDF_MION.resultKMeans(iK).SU_indCluster = IDX_SU;
    ClusteringSDF_MION.resultKMeans(iK).SU_sumD = SUMD_SU;
    ClusteringSDF_MION.resultKMeans(iK).time_indCluster = IDX_t;
    ClusteringSDF_MION.resultKMeans(iK).time_sumD = SUMD_t;
    
    % Cluster single units based on time series in TR resolution
    [IDX_SUTR, C, SUMD_SUTR] = kmeans(matFR_TR', K, 'Replicates', numReplicates, 'Options', opts);     
    
    ClusteringSDF_TR.resultKMeans(iK).SU_indCluster = IDX_SUTR;
    ClusteringSDF_TR.resultKMeans(iK).SU_sumD = SUMD_SUTR;
    
end

[a, c, totalSS] = kmeans(matNeuralRGR', 1);
ClusteringSDF_MION.totalSS = totalSS;
[a, c, totalSS] = kmeans(matFR_TR', 1);
ClusteringSDF_TR.totalSS = totalSS;

ClusteringSDF_MION.methods = 'KMeans';
ClusteringSDF_MION.setK = setK;
ClusteringSDF_MION.numReplicates = numReplicates;

ClusteringSDF_TR.methods = 'KMeans';
ClusteringSDF_TR.setK = setK;
ClusteringSDF_TR.numReplicates = numReplicates;

% (2) More fine temporal scale
% matFR_SU = NaN(375, length(validC));
FR_dTfine = createCellRegressor_indMov_discreteTime(dirDataNeural, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
    setMovie, 0.1); % in 10Hz (number of spikes for every 100ms)
% concatenate across movies
matFR_SU=[];
for iUnit = 1:size(FR_dTfine,1)
    tempFR = cat(1, FR_dTfine(iUnit, :).mnFR);
    matFR_SU(:,iUnit) = tempFR;
end
% (2-1) Normalized (z-scored)
matFR_SU_norm = zscore(matFR_SU); % noramlized time series

setK = 2:15;
opts = statset('Display','final');
numReplicates = 100; %5;

ClusteringSDF.methods = 'KMeans';
ClusteringSDF.setK = setK;
ClusteringSDF.numReplicates = numReplicates;

ClusteringSDFnorm.methods = 'KMeans';
ClusteringSDFnorm.setK = setK;
ClusteringSDFnorm.numReplicates = numReplicates;

% matCluster_SU = NaN(size(matR_SU, 2), 
for iK = 1:length(setK)
    
    K = setK(iK);
    % Cluster single units based on time series
    [IDX_SU, C, SUMD_SU] = kmeans(matFR_SU', K, 'Replicates', numReplicates, 'Options', opts); 
    % Cluster time points based on 50 singel unit correlation
    [IDX_t, C, SUMD_t] = kmeans(matFR_SU, K, 'Replicates', numReplicates, 'Options', opts);
    
    % Cluster single units based on normalized (z-scored) time series
    [IDX_SUnorm, C, SUMD_SUnorm] = kmeans(matFR_SU_norm', K, 'Replicates', numReplicates, 'Options', opts); 
    
    
    ClusteringSDF.resultKMeans(iK).SU_indCluster = IDX_SU;
    ClusteringSDF.resultKMeans(iK).SU_sumD = SUMD_SU;
    ClusteringSDF.resultKMeans(iK).time_indCluster = IDX_t;
    ClusteringSDF.resultKMeans(iK).time_sumD = SUMD_t;
    
    ClusteringSDFnorm.resultKMeans(iK).SU_indCluster = IDX_SUnorm;
    ClusteringSDFnorm.resultKMeans(iK).SU_sumD = SUMD_SUnorm;
    
end

[a, c, totalSS] = kmeans(matFR_SU', 1);
ClusteringSDF.totalSS = totalSS;
[a, c, totalSS] = kmeans(matFR_SU_norm', 1);
ClusteringSDFnorm.totalSS = totalSS;

% % plot Clustering results
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% for iK=1:length(setK)
%     plot(iK+1, ClusteringSDF.resultKMeans(iK).SU_sumD, 'ko', 'MarkerSize', 10, 'LineWidth', 2)
%     hold on
% end
% xlim([1 16])
% set(gca, 'LineWidth', 2, 'FontSize', 12)
% box off
% title('Clustering of single units based on time series')
% xlabel('Number of cluster')
% ylabel('Within-cluster distance')
% print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_SUSDF'), '-depsc');
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% for iK=1:length(setK)
%     plot(iK+1, ClusteringSDF.resultKMeans(iK).time_sumD, 'ko', 'MarkerSize', 10, 'LineWidth', 2)
%     hold on
% end
% xlim([1 16])
% set(gca, 'LineWidth', 2, 'FontSize', 12)
% box off
% title('Clustering of time points')
% xlabel('Number of cluster')
% ylabel('Within-cluster distance')
% print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_time'), '-depsc');


save(fullfile(dirDataNeural, sprintf('ClusteringSDF_%sMovie123.mat', nameSubjNeural)), 'ClusteringSDF*')


