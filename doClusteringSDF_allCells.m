function [] =  doClusteringSDF_allCells(flagParallel, flagSave)
% 2019/06/27 SHP
% Apply K-means clustering to time series of single units from AF and AM patches across multiple
% subjects


%% Settings
flagBiowulf = flagParallel;

if flagBiowulf
    dirDataHome = '/data/parks20/procdata/NeuroMRI/';
    %     addpath('/data/parks20/analysis/NeuroMRI/'); % to use doConv.m function
    dirDataNeural = fullfile(dirDataHome, 'Spi');
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/PROJECTS';
        dirProcdata = '/Volumes/PROCDATA';
        dirDataHome = fullfile(dirProcdata, 'parksh', '_macaque');
        dirLibrary = '/Volumes/LIBRARY';
        %         addpath(fullfile(dirLibrary, 'matlab_utils'));
        %         dirDataNeural = fullfile(dirDataHome, 'Spi');
    else % on virtual machine
        dirProjects = '/projects';
        dirProcdata = '/procdata';
        dirDataHome = fullfile(dirProcdata, 'parksh', '_macaque');
        dirLibrary = '/library';
        %         addpath(fullfile(dirLibrary, 'matlab_utils'));
        %         dirDataNeural = fullfile(dirDataHome, 'Spi');
    end
end

% Load the cocatenated time courses in different resolution
load(fullfile(dirDataHome, 'matSDF_Movie123_allCells.mat'), 'matSDF')


% Prepare the concatenate SDFs
matFR_TR = cat(2, matSDF.matFR_TR);
matNeuralRGR = cat(2, matSDF.matNeuralRGR);
matFR_SU = cat(2, matSDF.matFR_SU_10hz);
matFR_SU_norm = cat(2, matSDF.matFR_SU_10hz_norm);

% take care of the NaNs from matNeuralRGR
indValid = ~isnan(matNeuralRGR);
matNeuralRGR_noNaN = reshape(matNeuralRGR(indValid), 375-21, size(matNeuralRGR, 2));
clear matNeuralRGR
matNeuralRGR = matNeuralRGR_noNaN;
clear matNeuralRGR_noNaN


% Channel IDs for each subject
validChanIDs={}; validChanIDs_subj={};
[validChanIDs{1, 1:length(matSDF)}] = deal(matSDF.setCellIDs);
[validChanIDs_subj{1, 1:length(matSDF)}] = deal(matSDF.setCellIDs_subjID);
paramClustering_global.validChanIDs = validChanIDs;
paramClustering_global.validChanIDs_subjID = validChanIDs_subj;

%% K-means Clustering based on timeseries
% numRepeat = 2; %10; %100; % number of repetition for entire clustering

setK = 2:20; %15; %20; %15;
opts = statset('Display','final');
numReplicates = 100; %5; %

paramClustering_global.methods = 'KMeans';
paramClustering_global.setK = setK;
paramClustering_global.numReplicates = numReplicates;
paramClustering_global.descriptions = 'Clustering of neurons based on the masked whole-brain correlation maps: multiple repetitions (numRepeat) of an execution of kmeans function, which had its own multiple "replicates" (numReplicates)';
% paramClustering_global.numRepeat = numRepeat;

if flagParallel
    pool = parpool;                      % Invokes workers
    stream = RandStream('mlfg6331_64');  % Random number stream
    opts = statset('UseParallel',1,'UseSubstreams',1, 'Streams',stream,...
        'MaxIter', 1000, 'Display','final');
    paramClustering_global.parallel = opts;
end

% Total SS
[a, c, totalSS] = kmeans(matFR_TR', 1);
ClusteringSDF_TR.totalSS = totalSS;
[a, c, totalSS] = kmeans(matNeuralRGR', 1);
ClusteringSDF_MION.totalSS = totalSS;

[a, c, totalSS] = kmeans(matFR_SU', 1);
ClusteringSDF.totalSS = totalSS;
[a, c, totalSS] = kmeans(matFR_SU_norm', 1);
ClusteringSDFnorm.totalSS = totalSS;



for iK = 1:length(setK)
    
    K = setK(iK);
    
    
    fprintf(1, ':: K = %d; SDF ::\n', K);
    
    % 1. No MION
    [IDX_SUTR, C, SUMD_SUTR] = kmeans(matFR_TR', K, 'Replicates', numReplicates, 'Options', opts);
    [IDX_tTR, C, SUMD_tTR] = kmeans(matFR_TR, K, 'Replicates', numReplicates, 'Options', opts);
    ClusteringSDF_TR.resultKMeans(iK).SU_indCluster = IDX_SUTR;
    ClusteringSDF_TR.resultKMeans(iK).SU_sumD = SUMD_SUTR;
    ClusteringSDF_TR.resultKMeans(iK).time_indCluster = IDX_tTR;
    ClusteringSDF_TR.resultKMeans(iK).time_sumD = SUMD_tTR;
    
    % 2. MION
    [IDX_SUMION, C, SUMD_SUMION] = kmeans(matNeuralRGR', K, 'Replicates', numReplicates, 'Options', opts);
    [IDX_tMION, C, SUMD_tMION] = kmeans(matNeuralRGR, K, 'Replicates', numReplicates, 'Options', opts);
    ClusteringSDF_MION.resultKMeans(iK).SU_indCluster = IDX_SUMION;
    ClusteringSDF_MION.resultKMeans(iK).SU_sumD = SUMD_SUMION;
    ClusteringSDF_MION.resultKMeans(iK).time_indCluster = IDX_tMION;
    ClusteringSDF_MION.resultKMeans(iK).time_sumD = SUMD_tMION;
    
    % 3. Fine temporal resolution
    [IDX_SUFR, C, SUMD_SUFR] = kmeans(matFR_SU', K, 'Replicates', numReplicates, 'Options', opts);
    [IDX_tFR, C, SUMD_tFR] = kmeans(matFR_SU, K, 'Replicates', numReplicates, 'Options', opts);
    ClusteringSDF.resultKMeans(iK).SU_indCluster = IDX_SUFR;
    ClusteringSDF.resultKMeans(iK).SU_sumD = SUMD_SUFR;
    ClusteringSDF.resultKMeans(iK).time_indCluster = IDX_tFR;
    ClusteringSDF.resultKMeans(iK).time_sumD = SUMD_tFR;
    
    % 4. Fine temporal resolution: Normalized
    [IDX_SUFRnorm, C, SUMD_SUFRnorm] = kmeans(matFR_SU_norm', K, 'Replicates', numReplicates, 'Options', opts);
    [IDX_tFRnorm, C, SUMD_tFRnorm] = kmeans(matFR_SU_norm, K, 'Replicates', numReplicates, 'Options', opts);
    ClusteringSDFnorm.resultKMeans(iK).SU_indCluster = IDX_SUFRnorm;
    ClusteringSDFnorm.resultKMeans(iK).SU_sumD = SUMD_SUFRnorm;
    ClusteringSDFnorm.resultKMeans(iK).time_indCluster = IDX_tFRnorm;
    ClusteringSDFnorm.resultKMeans(iK).time_sumD = SUMD_tFRnorm;
    
    if flagSave
        save(fullfile(dirDataHome, 'ClusteringSDF_Movie123_allCells.mat'), ... %nameSubjNeural, nameSubjBOLD)),...
            'Clustering*', 'paramClustering*');
        fprintf(1, ':: K = %d; SDF clustering :: Results saved \n\n', K);
    end
end

% save(fullfile(dirDataNeural, sprintf('ClusteringSDF_%sMovie123.mat', nameSubjNeural)), 'ClusteringSDF*', 'paramClustering*')


%
% % matCluster_SU = NaN(size(matR_SU, 2),
% for iK = 1:length(setK)
%
%     K = setK(iK);
%
%     ClusteringSDF_MION.resultKMeans(iK).time_indCluster = IDX_t;
%     ClusteringSDF_MION.resultKMeans(iK).time_sumD = SUMD_t;
%
%
%
% end
%
% [a, c, totalSS] = kmeans(matNeuralRGR', 1);
% ClusteringSDF_MION.totalSS = totalSS;
% [a, c, totalSS] = kmeans(matFR_TR', 1);
% ClusteringSDF_TR.totalSS = totalSS;
%
% ClusteringSDF_MION.methods = 'KMeans';
% ClusteringSDF_MION.setK = setK;
% ClusteringSDF_MION.numReplicates = numReplicates;
%
% ClusteringSDF_TR.methods = 'KMeans';
% ClusteringSDF_TR.setK = setK;
% ClusteringSDF_TR.numReplicates = numReplicates;
%
% % (2) More fine temporal scale
% % matFR_SU = NaN(375, length(validC));
% FR_dTfine = createCellRegressor_indMov_discreteTime(dirDataNeural, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
%     setMovie, 0.1); % in 10Hz (number of spikes for every 100ms)
% % concatenate across movies
% matFR_SU=[];
% for iUnit = 1:size(FR_dTfine,1)
%     tempFR = cat(1, FR_dTfine(iUnit, :).mnFR);
%     matFR_SU(:,iUnit) = tempFR;
% end
% % (2-1) Normalized (z-scored)
% matFR_SU_norm = zscore(matFR_SU); % noramlized time series
%
% setK = 2:15;
% opts = statset('Display','final');
% numReplicates = 100; %5;
%
% ClusteringSDF.methods = 'KMeans';
% ClusteringSDF.setK = setK;
% ClusteringSDF.numReplicates = numReplicates;
%
% ClusteringSDFnorm.methods = 'KMeans';
% ClusteringSDFnorm.setK = setK;
% ClusteringSDFnorm.numReplicates = numReplicates;
%
% % matCluster_SU = NaN(size(matR_SU, 2),
% for iK = 1:length(setK)
%
%     K = setK(iK);
%     % Cluster single units based on time series
%     [IDX_SUMION, C, SUMD_SUMION] = kmeans(matFR_SU', K, 'Replicates', numReplicates, 'Options', opts);
%     % Cluster time points based on 50 singel unit correlation
%     [IDX_t, C, SUMD_t] = kmeans(matFR_SU, K, 'Replicates', numReplicates, 'Options', opts);
%
%     % Cluster single units based on normalized (z-scored) time series
%     [IDX_SUFRnorm, C, SUMD_SUnorm] = kmeans(matFR_SU_norm', K, 'Replicates', numReplicates, 'Options', opts);
%
%
%     ClusteringSDF.resultKMeans(iK).SU_indCluster = IDX_SUMION;
%     ClusteringSDF.resultKMeans(iK).SU_sumD = SUMD_SUMION;
%     ClusteringSDF.resultKMeans(iK).time_indCluster = IDX_t;
%     ClusteringSDF.resultKMeans(iK).time_sumD = SUMD_t;
%
%     ClusteringSDFnorm.resultKMeans(iK).SU_indCluster = IDX_SUFRnorm;
%     ClusteringSDFnorm.resultKMeans(iK).SU_sumD = SUMD_SUnorm;
%
% end
%
% [a, c, totalSS] = kmeans(matFR_SU', 1);
% ClusteringSDF.totalSS = totalSS;
% [a, c, totalSS] = kmeans(matFR_SU_norm', 1);
% ClusteringSDFnorm.totalSS = totalSS;

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



