function [] =  doClusteringSDF_multipleSubj_prob(flagParallel, flagSave)
% 2017/02/02 SHP
% Apply K-means clustering to time series of singel units from multiple
% subjects
% Repeat the k-means procedure multiple times to get likelihood

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
        dirDataHome = fullfile(dirProcdata, 'parksh');
        dirLibrary = '/Volumes/LIBRARY';
%         addpath(fullfile(dirLibrary, 'matlab_utils'));
        dirDataNeural = fullfile(dirDataHome, 'Spi');
    else % on virtual machine
        dirProjects = '/projects';
        dirProcdata = '/procdata';
        dirDataHome = fullfile(dirProcdata, 'parksh');
        dirLibrary = '/library';
%         addpath(fullfile(dirLibrary, 'matlab_utils'));
        dirDataNeural = fullfile(dirDataHome, 'Spi');
    end
end

%% Prepare time-series matrices in different resolution
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
setMovie = [1 2 3];
% nt = 375;
% 
% % MION function
% TR=2.4;
% k = gampdf([-40:TR:40],4,2);
% 
% 
% matSDF = struct([]);
% for iSubj = 1:length(setNameSubjNeural)
%     
%     % Load the data
%     nameSubjNeural = setNameSubjNeural{iSubj};
%     dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
%     filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
%     load(fullfile(dirDataNeural, filenameNeural))
%     fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
%         
%     switch lower(nameSubjNeural)
%         case 'spi'
%             excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
%             validC = setdiff(1:length(paramSDF.setCellIDs), excChanIndex)';
%         otherwise
%             [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
%             validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
%     end   
%     matSDF(iSubj).setCellIDs = paramSDF.setCellIDs(validC);
%     
%     indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
%     % 1. No MION, TR resolution (375 time points - 7 points at the )
%     matFR_TR = NaN(nt, length(validC));
%     for iChan = 1:length(validC) 
%         tempFR = [];
%         tempFR = cat(1, S(validC(iChan), indMovieNeuron).mnFR);
%         matFR_TR(:,iChan) = tempFR;
%     end
%     matSDF(iSubj).matFR_TR = matFR_TR;
%     
%     % 2. MION convolved, TR resolution
%     matNeuralRGR = NaN(nt, length(validC));
%     for iChan = 1:length(validC) 
%         neuralrgrs=[];
%         neuralrgrs = cat(1, S(validC(iChan), indMovieNeuron).mnFR);
%         neuralrgrs = neuralrgrs-mean(neuralrgrs); % centering
%         neuralrgrs = doConv(neuralrgrs,k); % convolve MI
%         matNeuralRGR(:,iChan) = neuralrgrs';
%     end
%     matSDF(iSubj).matNeuralRGR = matNeuralRGR;
%     
%     % 3. Fine temporal resolution
%     FR_dTfine = createCellRegressor_indMov_discreteTime(dirDataNeural, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
%         setMovie, 0.1); % in 10Hz (number of spikes for every 100ms)
%     % concatenate across movies
%     matFR_SU=[];
%     for iUnit = 1:size(FR_dTfine,1)
%         tempFR = cat(1, FR_dTfine(iUnit, :).mnFR);
%         matFR_SU(:,iUnit) = tempFR;
%     end
%     matSDF(iSubj).matFR_SU = matFR_SU;
%     
%     % 4. Fine temporal resolution, normalized (z-scored)
%     matFR_SU_norm = zscore(matFR_SU); % noramlized time series
%     matSDF(iSubj).matFR_SU_norm = matFR_SU_norm;
%     
% end

% load the cocatenated time courses in different resolution
load(fullfile(dirDataNeural, sprintf('matSDF_%sMovie123.mat', cell2mat(setNameSubjNeural))), 'matSDF')


% Concatenate SDFs 
matFR_TR = cat(2, matSDF.matFR_TR);
matNeuralRGR = cat(2, matSDF.matNeuralRGR);
matFR_SU = cat(2, matSDF.matFR_SU);
matFR_SU_norm = cat(2, matSDF.matFR_SU_norm);

% Channel IDs for each subject
validChanIDs_subj={};
[validChanIDs_subj{1, 1:length(matSDF)}] = deal(matSDF.setCellIDs);
paramClustering_global.validChanIDs_Subj = validChanIDs_subj;

%% K-means Clustering based on timeseries 
numRepeat = 2; %10; %100; % number of repetition for entire clustering

setK = 2:15; %20; %15;
opts = statset('Display','final');
numReplicates = 5; %

paramClustering_global.methods = 'KMeans';
paramClustering_global.setK = setK;
paramClustering_global.numReplicates = numReplicates;
paramClustering_global.descriptions = 'Clustering of neurons based on the masked whole-brain correlation maps: multiple repetitions (numRepeat) of an execution of kmeans function, which had its own multiple "replicates" (numReplicates)';
paramClustering_global.numRepeat = numRepeat;


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
%     
%     % 1. No MION
%     SU_indCluster_FR_TR = NaN(size(matFR_TR, 2), numRepeat);
%     SU_sumD_FR_TR = NaN(K, numRepeat);
%     
%     % 2. MION
%     SU_indCluster_NeuralRGR = NaN(size(matNeuralRGR, 2), numRepeat);
%     SU_sumD_NeuralRGR = NaN(K, numRepeat);
    
    % 3. Fine temporal resolution
    SU_indCluster_FR_SU = NaN(size(matFR_SU, 2), numRepeat);
    SU_sumD_FR_SU = NaN(K, numRepeat);
    
    % 3. Fine temporal resolution: Normalized
    SU_indCluster_FR_SU_norm = NaN(size(matFR_SU_norm, 2), numRepeat);
    SU_sumD_FR_SU_norm = NaN(K, numRepeat);
    

    for iRep = 1:numRepeat
        
        fprintf(1, ':: K = %d; SDF ::\n', K);
                       
%         % 1. No MION
%         [IDX_SUTR, C, SUMD_SUTR] = kmeans(matFR_TR', K, 'Replicates', numReplicates, 'Options', opts);
%         
%         SU_indCluster_FR_TR(:, iRep) = IDX_SUTR;
%         SU_sumD_FR_TR(:, iRep) = SUMD_SUTR;
%         
%         % 2. MION
%         [IDX_SUMION, C, SUMD_SUMION] = kmeans(matNeuralRGR', K, 'Replicates', numReplicates, 'Options', opts);
%         
%         SU_indCluster_NeuralRGR(:, iRep) = IDX_SUMION;
%         SU_sumD_NeuralRGR(:, iRep) = SUMD_SUMION;
        
        % 3. Fine temporal resolution
        [IDX_SUFR, C, SUMD_SUFR] = kmeans(matFR_SU', K, 'Replicates', numReplicates, 'Options', opts);
        
        SU_indCluster_FR_SU(:, iRep) = IDX_SUFR;
        SU_sumD_FR_SU(:, iRep) = SUMD_SUFR;
        
        % 4. Fine temporal resolution: Normalized
        [IDX_SUFRnorm, C, SUMD_SUFRnorm] = kmeans(matFR_SU_norm', K, 'Replicates', numReplicates, 'Options', opts);
        
        SU_indCluster_FR_SU_norm(:, iRep) = IDX_SUFRnorm;
        SU_sumD_FR_SU_norm(:, iRep) = SUMD_SUFRnorm;
   
    end

%     ClusteringSDF_TR.resultKMeans(iK).SU_indCluster = SU_indCluster_FR_TR;
%     ClusteringSDF_TR.resultKMeans(iK).SU_sumD = SU_sumD_FR_TR;
%     
%     ClusteringSDF_MION.resultKMeans(iK).SU_indCluster = SU_indCluster_NeuralRGR;
%     ClusteringSDF_MION.resultKMeans(iK).SU_sumD = SU_sumD_NeuralRGR;
    
    ClusteringSDF.resultKMeans(iK).SU_indCluster = SU_indCluster_FR_SU;
    ClusteringSDF.resultKMeans(iK).SU_sumD = SU_sumD_FR_SU;
    
    ClusteringSDFnorm.resultKMeans(iK).SU_indCluster = SU_indCluster_FR_SU_norm;
    ClusteringSDFnorm.resultKMeans(iK).SU_sumD = SU_sumD_FR_SU_norm;
    
    
%     if flagSave
%         save(fullfile(dirDataNeural, sprintf('ClusteringSDF_%sMovie123.mat', cell2mat(setNameSubjNeural))), ... %nameSubjNeural, nameSubjBOLD)),...
%             'Clustering*', 'paramClustering*');
%         fprintf(1, ':: K = %d; SDF clustering :: Results saved \n\n', K);
%     end
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



