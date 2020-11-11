function [] =  doClusteringSDF_allCells(flagParallel, flagSave)
%
% 2020/10/05 SHP
% Modify the code to match to the data structure from the new version of
% "saveMovieData_SU_SDF_allCells.m"
% 2019/06/27 SHP
% Apply K-means clustering to time series of single units from AF and AM patches across multiple
% subjects


% %%%%%%%% code in progress %%%
% Normalization: z-score? using separate baseline? before clustering?
% K-means distance metric: squared Euclidean? or correlation?
% order can be done using within-cluster distance, or distance matrix with
% clustering


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
load(fullfile(dirDataHome, 'matSDF_Movie123_allCells.mat'), 'matTS_FP')


% Prepare the concatenate SDFs (time x cell)
matFR_TR = cat(2, matTS_FP.matFR_TR);
matNeuralRGR = cat(2, matTS_FP.matNeuralRGR);
matFR_SU_10hz = cat(2, matTS_FP.matFR_SU_10hz);
matFR_SU_1hz = cat(2, matTS_FP.matFR_SU_1hz);

% take care of the NaNs from matNeuralRGR
indValid = ~isnan(matNeuralRGR);
matNeuralRGR_noNaN = reshape(matNeuralRGR(indValid), 375-21, size(matNeuralRGR, 2));
clear matNeuralRGR
matNeuralRGR = matNeuralRGR_noNaN;
clear matNeuralRGR_noNaN


% % Channel IDs for each subject
% validChanIDs={}; validChanIDs_subj={};
% [validChanIDs{1, 1:length(matSDF)}] = deal(matSDF.setCellIDs);
% [validChanIDs_subj{1, 1:length(matSDF)}] = deal(matSDF.setCellIDs_subjID);
% paramClustering_global.validChanIDs = validChanIDs;
% paramClustering_global.validChanIDs_subjID = validChanIDs_subj;

%% K-means Clustering based on timeseries
% numRepeat = 2; %10; %100; % number of repetition for entire clustering

setK = 2:20; %15; %20; %15;
opts = statset('Display','final');
numReplicates = 5; %100; %5; %

paramClustering_global.methods = 'KMeans';
paramClustering_global.setK = setK;
paramClustering_global.numReplicates = numReplicates;
paramClustering_global.descriptions = 'Clustering of neurons based on the time series multiple repetitions (numRepeat) of an execution of kmeans function, which had its own multiple "replicates" (numReplicates)';
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
% [a, c, totalSS] = kmeans(matNeuralRGR', 1);
% ClusteringSDF_MION.totalSS = totalSS;

[a, c, totalSS] = kmeans(matFR_SU_10hz', 1);
ClusteringSDF_10hz.totalSS = totalSS;
[a, c, totalSS] = kmeans(matFR_SU_1hz', 1);
ClusteringSDF_1hz.totalSS = totalSS;



for iK = 1:length(setK)
    
    K = setK(iK);
    
    
    fprintf(1, ':: K = %d; SDF ::\n', K);
    
%     % 1. No MION
%     [IDX_SUTR, C, SUMD_SUTR] = kmeans(matFR_TR', K, 'Replicates', numReplicates, 'Options', opts);
%     [IDX_tTR, C, SUMD_tTR] = kmeans(matFR_TR, K, 'Replicates', numReplicates, 'Options', opts);
%     ClusteringSDF_TR.resultKMeans(iK).SU_indCluster = IDX_SUTR;
%     ClusteringSDF_TR.resultKMeans(iK).SU_sumD = SUMD_SUTR;
%     ClusteringSDF_TR.resultKMeans(iK).time_indCluster = IDX_tTR;
%     ClusteringSDF_TR.resultKMeans(iK).time_sumD = SUMD_tTR;
%     
%     % 2. MION
%     [IDX_SUMION, C, SUMD_SUMION] = kmeans(matNeuralRGR', K, 'Replicates', numReplicates, 'Options', opts);
%     [IDX_tMION, C, SUMD_tMION] = kmeans(matNeuralRGR, K, 'Replicates', numReplicates, 'Options', opts);
%     ClusteringSDF_MION.resultKMeans(iK).SU_indCluster = IDX_SUMION;
%     ClusteringSDF_MION.resultKMeans(iK).SU_sumD = SUMD_SUMION;
%     ClusteringSDF_MION.resultKMeans(iK).time_indCluster = IDX_tMION;
%     ClusteringSDF_MION.resultKMeans(iK).time_sumD = SUMD_tMION;
%     
%     % 3. Fine temporal resolution: 10 hz
%     [IDX_SUFR, C, SUMD_SUFR] = kmeans(matFR_SU_10hz', K, 'Replicates', numReplicates, 'Options', opts);
%     [IDX_tFR, C, SUMD_tFR] = kmeans(matFR_SU_10hz, K, 'Replicates', numReplicates, 'Options', opts);
%     ClusteringSDF_10hz.resultKMeans(iK).SU_indCluster = IDX_SUFR;
%     ClusteringSDF_10hz.resultKMeans(iK).SU_sumD = SUMD_SUFR;
%     ClusteringSDF_10hz.resultKMeans(iK).time_indCluster = IDX_tFR;
%     ClusteringSDF_10hz.resultKMeans(iK).time_sumD = SUMD_tFR;
    
    % 4. Fine temporal resolution: 1hz
    [IDX_SUFR_1hz, C, SUMD_SUFR_1hz] = kmeans(matFR_SU_1hz', K, 'Replicates', numReplicates, 'Options', opts, 'Distance', 'correlation');
    [IDX_tFR_1hz, C, SUMD_tFR_1hz] = kmeans(matFR_SU_1hz, K, 'Replicates', numReplicates, 'Options', opts, 'Distance', 'correlation');
    ClusteringSDF_1hz.resultKMeans(iK).SU_indCluster = IDX_SUFR_1hz;
    ClusteringSDF_1hz.resultKMeans(iK).SU_sumD = SUMD_SUFR_1hz;
    ClusteringSDF_1hz.resultKMeans(iK).time_indCluster = IDX_tFR_1hz;
    ClusteringSDF_1hz.resultKMeans(iK).time_sumD = SUMD_tFR_1hz;
    
    if flagSave
        save(fullfile(dirDataHome, 'ClusteringSDF_Movie123_allCells.mat'), ... %nameSubjNeural, nameSubjBOLD)),...
            'Clustering*', 'paramClustering*');
        fprintf(1, ':: K = %d; SDF clustering :: Results saved \n\n', K);
    end
end



%% quick check
setK = paramClustering_global.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
%     matWSS_corrMap(:,iK) = sum(Clustering_moviemask_valid.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
    
    matWSS(:,iK) = sum(ClusteringSDF_1hz.resultKMeans(iK).SU_sumD); %
%     matWSS_SDFnorm(:,iK) = sum(ClusteringSDFnorm.resultKMeans(iK).SU_sumD); %
%     matWSS_SDF_TR(:,iK) = sum(ClusteringSDF_TR.resultKMeans(iK).SU_sumD); %
%     matWSS_SDF_MION(:,iK) = sum(ClusteringSDF_MION.resultKMeans(iK).SU_sumD); %
end

totalSS= ClusteringSDF_1hz.totalSS;
% betweenSS_corrMap = totalSS_corrMap-matWSS_corrMap;
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

figure;
plot(setK, propExplained, 'o-')

iK = 15; %19;
[sortedClust, indSortChan] = sort(ClusteringSDF_1hz.resultKMeans(iK).SU_indCluster);

cMap_Area = [91 148 203; 237 28 35; 248 148 29; 6 177 102]./255; % from Kenji's schematic
cMap_Area(4, :) = cMap_Area(4, :).*0.7; % make the green a bit darker

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 400 910 675])

sp1 = subplot('Position', [0.2 0.2 0.7 0.7]);
imagesc(zscore(matTS_FP.matFR_SU_1hz(:, indSortChan))); %(indSortChan, orderROI)')
set(sp1, 'CLim', [-1 1].*5)
% set(sp1, 'YTick', 1:numROI, 'YTickLabel', Clustering_meanROI.nameROI(orderROI))
locDiff = cat(1, 0, find(diff(sortedClust)>0), length(sortedClust));
set(sp1, 'XTick', locDiff+0.5, 'XTickLabel', locDiff)
title(sprintf('Clustered cells from 4 FPs using 1 hz time courses: K=%d', iK+1))
% xlabel('Cumulative number of cells')
colormap(sp1, hot)
set(gca, 'TickDir', 'out', 'Box', 'off')
line([locDiff+0.5 locDiff+0.5]', repmat(get(gca, 'YLim')', 1, length(locDiff)), 'Color', 'k', 'LineWidth', 0.5)

set(sp1, 'XTickLabel', [])
spp1 = axes('Position', sp1.Position, 'Color', 'none', 'XAxisLocation', 'top');
spp1.XLim = sp1.XLim;
set(spp1, 'XTick', locDiff+0.5, 'XTickLabel', locDiff, 'TickDir', 'out', 'XTickLabel', [], 'Box', 'off', 'YColor', 'none')

sp2 = subplot('Position', [0.2 0.08 0.7 0.05]);
imagesc(matTS_FP.catAreaID(indSortChan)') %(indSortChan)')
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



