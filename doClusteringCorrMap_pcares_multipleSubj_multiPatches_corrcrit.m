function [] = doClusteringCorrMap_pcares_multipleSubj_multiPatches_corrcrit(flagParallel, flagSave)
% 2019/07/08 SHP
% Apply K-means clustering to pca-res correlation map of single units from multiple
% subjects from multiple patches, using only a subset of voxels that
% survive a criterion of correlation
% 


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
    
% Add necessary toolbox 
addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution

% Set directories 
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was'};
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh/_macaque');
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% % load the masks (movie-driven & brain-only mask)
% load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');

% setArea = {'AF', 'AM', 'AM+', 'ML'};
% setNameSubjNeural{1} = {'Tor', 'Rho', 'Sig', 'Spi', 'Moc'}; % for AF: 11 12 13 14 15
% setNameSubjNeural{2} = {'Mat', 'Was'}; % for AM: 21 22
% setNameSubjNeural{3} = {'Dan', 'Moc'}; % for AM+: 31 32

%% Load the corr map and select valid channel: subject by subject
numSubject = size(setNameSubjNeural, 2);
matR_SU_all = [];

for iSubj = 1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_pcares_masked.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');
    
    % For each cell, save the subject ID: Mochi's AF cells get 15, Mochi's
    % AM cells get 32
    switch lower(nameSubjNeural) % for AF: 11 12 13 14 15; for AM: 21 22; for AM+: 31 32
        case 'tor'
            subjID = 11;
        case 'rho'
            subjID = 12;
        case 'sig'
            subjID = 13;
        case 'spi'
            subjID = 14;
        case 'mat' %AM
            subjID = 21;
        case 'dan' %AM+
            subjID = 31;
        case 'moc' % AF and AM+
            subjID = [15 32];
        case 'was' %AM
            subjID = 22;
    end
        
    paramClustering(iSubj).nameSubj = nameSubjNeural;
    paramClustering(iSubj).validChanIndex = (1:length(paramCorr.validChanIndex))'; %validChanIndex_clustering;
    paramClustering(iSubj).validChanID = paramCorr.validChanID; %validChanID_clustering;
    paramClustering(iSubj).validChan_subjID = ones(size(paramClustering(iSubj).validChanIndex)).*subjID(1);
     if strcmpi(nameSubjNeural, 'moc') % Mochi has both AF and AM+ cells, so assign different IDs for cells from each area
        locAM = ~cellfun(@isempty, strfind(cellstr(paramCorr.validChanID), 'AM'));
        paramClustering(iSubj).validChan_subjID(locAM) = subjID(2);
     end
    
    matR_SU_valid = matR_SU(:, paramClustering(iSubj).validChanIndex);
    matR_SU_all = cat(2, matR_SU_all, matR_SU_valid);
    
     clear matR_SU matR_SU_valid
     
         
%     % 2017/01/19: Excluded some units, in case of Spice's data which hasn't been went through this procedure yet
%     switch lower(nameSubjNeural)
%         case 'spi'
%             % 2018 Jan dataset
%             excChanIndex = [6 11 30 35]; % for 2018 Jan dataset: cells that are not face selective
%             validChanIndex_clustering = setdiff(1:length(paramCorr.validChanIndex), excChanIndex)';
%             validChanID_clustering = cat(2, paramCorr.validChanID(validChanIndex_clustering,:), num2str(ones(size(validChanIndex_clustering)).*iSubj)); %paramCorr.validChanID(validChanIndex_clustering,:);
% 
%             %           % 2016 Nov Dataset
% %             excChanIndex = [10 13 22 27 30 49]; % 2016 Nov dataset: cells were not same acrossd two days
% %             validChanIndex_clustering = setdiff(paramCorr.validChanIndex, excChanIndex);
% %             validChanID_clustering = cat(2, paramCorr.validChanID(validChanIndex_clustering,:), num2str(ones(size(validChanIndex_clustering)).*iSubj)); %paramCorr.validChanID(validChanIndex_clustering,:);
%         case 'dav'
%             validChanIndex_clustering = [1 4 5 8 9 10]'; % setdiff(paramCorr.validChanIndex, excChanIndex);
%             validChanID_clustering = cat(2, paramCorr.validChanID(validChanIndex_clustering,:), num2str(ones(size(validChanIndex_clustering)).*iSubj));
%         case 'dan'
%             validChanIndex_clustering = [3 4 5 6 7 8 9]'; % setdiff(paramCorr.validChanIndex, excChanIndex);
%             validChanID_clustering = cat(2, paramCorr.validChanID(validChanIndex_clustering,:), num2str(ones(size(validChanIndex_clustering)).*iSubj));
%         otherwise
%             validChanIndex_clustering = (1:length(paramCorr.validChanIndex))';
%             validChanID_clustering = cat(2, paramCorr.validChanID, num2str(ones(size(paramCorr.validChanIndex)).*iSubj)); %paramCorr.validChanID;
%     end
end


%% Clustering on masked map 
% [nx ny nz] = size(movieDrivenAmp.mask_amp1);
% nVox = nx*ny*nz;

% PCA-residual maps are already movie-masked
% moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D
matR_SU_all_moviemask = matR_SU_all; %matR_SU_all(moviemask_vec,:); %matR_SU(moviemask_vec,:); % 15495 voxels
% brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
% matR_SU_all_brainmask = matR_SU_all(brainmask_vec,:); %matR_SU(brainmask_vec,:); % 27113 voxels


%% Voxel Selection
%% critCorr1: voxels that showed pca-res correlation higher than 0.3 with any 5% of the entire population
% critCorr = 0.3; 
% critNumNeuron = ceil(size(matR_SU_all_moviemask, 2).*0.05); % 5% of the population %6; %13;
% matValidVox = abs(matR_SU_all_moviemask)>critCorr;
% locValidVox = sum(matValidVox, 2)>critNumNeuron;
% matR_SU_all_moviemask_valid  = matR_SU_all_moviemask(locValidVox, :);
% 
% paramClustering_global.critCorr = critCorr;
% paramClustering_global.critNumNeuron = critNumNeuron;
% paramClustering_global.locValidVox = locValidVox;

%% critCorr2: voxels that showed original correlation (not pca-res correlation) higher than 0.3 with any 5% of the entire population
%% For critCorr2, we use the same set of voxels as in "Clustering_%%sMovie123_new_masked_critCorr1.mat"
% load the parameters that used previously
load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)),...
    'paramClustering_global') 
locValidVox = paramClustering_global.locValidVox; 
paramClustering_global.locValidVox_description = 'same set of voxels that are used in original (not pcares) correlation map critCorr1';
matR_SU_all_moviemask_valid  = matR_SU_all_moviemask(locValidVox, :);

%% Clustering params
% numRepeat = 100; % number of repetition for entire clustering
setK = 2:20; %15;
opts = statset('Display','final');
numReplicates = 100; %15; %100; % 5; %10; %5;

paramClustering_global.methods = 'KMeans';
paramClustering_global.setK = setK;
paramClustering_global.numReplicates = numReplicates;
paramClustering_global.descriptions = 'Clustering of neurons based on the masked whole-brain correlations between each neuron and 1PC-regressed out fMRI responses';
% paramClustering_global.numRepeat = numRepeat;

if flagParallel
    pool = parpool;                      % Invokes workers
    stream = RandStream('mlfg6331_64');  % Random number stream
    opts = statset('UseParallel',1,'UseSubstreams',1, 'Streams',stream,...
        'MaxIter', 1000, 'Display','final');
    paramClustering_global.parallel = opts;
end

% totalSS
% [a, c, totalSS] = kmeans(matR_SU_all_moviemask', 1); % single unit
% Clustering_moviemask.totalSS_SU = totalSS;
[a, c, totalSS] = kmeans(matR_SU_all_moviemask_valid, 1); % voxel
Clustering_moviemask_valid.totalSS_Vox = totalSS;
% [a, c, totalSS] = kmeans(matR_SU_all_brainmask', 1);
% Clustering_brainmask.totalSS = totalSS;

for iK = 1:length(setK)
    
    K = setK(iK);                

    fprintf(1, ':: K = %d; Movie-driven mask : only valid voxels ::\n', K);
    
    % Cluster single units based on whole brain correlation
    [IDX_SU, C_SU, SUMD_SU] = kmeans(matR_SU_all_moviemask_valid', K,...
        'Replicates', numReplicates, 'Options', opts);

    % Cluster voxels based on single unit correlation
    [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_all_moviemask_valid, K,...
        'Replicates', numReplicates, 'Options', opts);

    
    Clustering_moviemask_valid.resultKMeans(iK).SU_indCluster = IDX_SU; %SU_indCluster_moviemask;
    Clustering_moviemask_valid.resultKMeans(iK).SU_sumD = SUMD_SU; %SU_sumD_moviemask;
    Clustering_moviemask_valid.resultKMeans(iK).Vox_indCluster = IDX_Vox; %Vox_indCluster_moviemask;
    Clustering_moviemask_valid.resultKMeans(iK).Vox_sumD = SUMD_Vox; %Vox_sumD_moviemask;
            
    if flagSave
        save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_pcares_masked_critCorr2.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
            'Clustering*', 'paramClustering*');
%         save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_pcares_masked.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
%             'Clustering*', 'paramClustering*');
        fprintf(1, ':: K = %d; Movie-driven mask : only valid voxels :: Results saved \n', K);
    end

end

%%%   
% plot Clustering results
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iK=1:length(setK)
    plot(iK+1, Clustering_moviemask_valid.resultKMeans(iK).SU_sumD, 'ko', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
end
xlim([1 16])
set(gca, 'LineWidth', 2, 'FontSize', 12)
box off
title('Clustering of single units')
xlabel('Number of cluster')
ylabel('Within-cluster distance')
% print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_SU'), '-depsc');

%% Perform gap statistics to get the optimal number of clusters
% Gap statistics to get the optimal number
fprintf(1, 'Performing gap statistics on SU clustering... \n');
eva_clusteringSU = evalclusters(matR_SU_all_moviemask_valid', 'kmeans', 'gap', 'KList', paramClustering_global.setK, 'Distance', 'sqEuclidean');
save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_pcares_masked_critCorr2_eval.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
        'eva_clusteringSU');
fprintf(1, '...done! Evaluation results saved. \n')

fprintf(1, 'Performing gap statistics on Voxel clustering... \n');
eva_clusteringVox = evalclusters(matR_SU_all_moviemask_valid, 'kmeans', 'gap', 'KList', paramClustering_global.setK, 'Distance', 'sqEuclidean'); 
fprintf(1, '...done! Evaluation results saved. \n')
save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_pcares_masked_critCorr2_eval.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
        'eva_clusteringSU', 'eva_clusteringVox');

if flagParallel
    delete(pool);
end

% % 2017/01/19 Record the results of each replicate (number of iterations,
% % sum of the distance)
% % diary(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked_ReplicatesRecords', nameSubjNeural, nameSubjBOLD)))
% for iK = 1:length(setK)
%     
%     K = setK(iK);   
%     
% %     SU_indCluster_moviemask = NaN(sum(locValidVox), numRepeat);
% %     SU_sumD_moviemask = NaN(K, numRepeat);
% %     
% %     Vox_indCluster_moviemask = NaN(sum(locValidVox), numRepeat);
% %     Vox_sumD_moviemask = NaN(K, numRepeat);
% % 
% %     for iRep = 1:numRepeat
%     
%         fprintf(1, ':: K = %d; Movie-driven mask ::\n', K);
%         
%         % Cluster single units based on whole brain correlation
%         [IDX_SU, C_SU, SUMD_SU] = kmeans(matR_SU_all_moviemask_valid', K,...
%             'Replicates', numReplicates, 'Options', opts);
%         
% %         SU_indCluster_moviemask(:, iRep) = IDX_SU;
% %         SU_sumD_moviemask(:, iRep) = SUMD_SU;
%         
%         % Cluster voxels based on 50 singel unit correlation
%         [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_all_moviemask_valid, K,...
%             'Replicates', numReplicates, 'Options', opts);
%         
% %         Vox_indCluster_moviemask(:, iRep) = IDX_Vox;
% %         Vox_sumD_moviemask(:, iRep) = SUMD_Vox;
%         
% %     end
% 
%     Clustering_moviemask_valid.resultKMeans(iK).SU_indCluster = IDX_SU; %SU_indCluster_moviemask;
%     Clustering_moviemask_valid.resultKMeans(iK).SU_sumD = SUMD_SU; %SU_sumD_moviemask;
%     Clustering_moviemask_valid.resultKMeans(iK).Vox_indCluster = IDX_Vox; %Vox_indCluster_moviemask;
%     Clustering_moviemask_valid.resultKMeans(iK).Vox_sumD = SUMD_Vox; %Vox_sumD_moviemask;
%     
%     if flagSave
% %         save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_pcares_masked_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
% %             'Clustering*', 'paramClustering*');
%         fprintf(1, ':: K = %d; Movie-driven mask :: Results saved \n', K);
%     end
% % 
% 
% %     Clustering_moviemask.resultKMeans(iK).SU_indCluster = IDX_SU;
% %     Clustering_moviemask.resultKMeans(iK).SU_sumD = SUMD_SU;
% %     Clustering_moviemask.resultKMeans(iK).SU_totalD = totalD;
% %     Clustering_moviemask.resultKMeans(iK).SU_centeroid = C_SU;
% %     Clustering_moviemask.resultKMeans(iK).Vox_indCluster = IDX_Vox;
% %     Clustering_moviemask.resultKMeans(iK).Vox_sumD = SUMD_Vox;
% 
% 
% end

% %% Perform gap statistics to get the optimal number of clusters
% % Gap statistics to get the optimal number
% fprintf(1, 'Performing gap statistics on SU clustering... \n');
% eva_clusteringSU = evalclusters(matR_SU_all_moviemask_valid', 'kmeans', 'gap', 'KList', paramClustering_global.setK, 'Distance', 'sqEuclidean');
% save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_pcares_masked_critCorr1_eval.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
%         'eva_clusteringSU');
% fprintf(1, '...done! Evaluation results saved. \n')
% 
% fprintf(1, 'Performing gap statistics on Voxel clustering... \n');
% eva_clusteringVox = evalclusters(matR_SU_all_moviemask_valid, 'kmeans', 'gap', 'KList', paramClustering_global.setK, 'Distance', 'sqEuclidean'); 
% fprintf(1, '...done! Evaluation results saved. \n')
% save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_pcares_masked_critCorr1_eval.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
%         'eva_clusteringSU', 'eva_clusteringVox');
% 
% if flagParallel
%     delete(pool);
% end

% for iK = 1:length(setK)
%     
%     K = setK(iK);             
%     
%     diary on
%     fprintf(1, ':: K = %d; Brain mask ::\n', K);
%     % 2. Brain mask
%     % Cluster single units based on whole brain correlation
%     [IDX_SU, C_SU, SUMD_SU, totalD] = kmeans(matR_SU_all_brainmask', K, 'Replicates', numReplicates, 'Options', opts); 
%     
%     diary off
%     % Cluster voxels based on 50 singel unit correlation
%     [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_all_brainmask, K, 'Replicates', 5, 'Options', opts); %, numReplicates, 'Options', opts);
%     
%     Clustering_brainmask.resultKMeans(iK).SU_indCluster = IDX_SU;
%     Clustering_brainmask.resultKMeans(iK).SU_sumD = SUMD_SU;
%     Clustering_brainmask.resultKMeans(iK).SU_totalD = totalD;
%     Clustering_brainmask.resultKMeans(iK).SU_centeroid = C_SU;
%     Clustering_brainmask.resultKMeans(iK).Vox_indCluster = IDX_Vox;
%     Clustering_brainmask.resultKMeans(iK).Vox_sumD = SUMD_Vox;
%     
% end
% save(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
%     'Clustering_moviemask', 'Clustering_brainmask', 'paramClustering*');

% %%%   
% % plot Clustering results
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% for iK=1:length(setK)
%     plot(iK+1, mean(Clustering_moviemask_valid.resultKMeans(iK).SU_sumD), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
%     hold on
% end
% xlim([1 16])
% set(gca, 'LineWidth', 2, 'FontSize', 12)
% box off
% title('Clustering of single units')
% xlabel('Number of cluster')
% ylabel('Within-cluster distance')
% % print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_SU'), '-depsc');


% %% Gap statistics
% tempRand = -1 + (1+1).*rand(size(matR_SU)); % uniformly distributed random samples in interval [-1, 1]
% for iK = 1:length(setK)
%     
%     K = setK(iK);
%     % Cluster single units based on whole brain correlation
%     [IDX_SU_Rand, C, SUMD_SU_Rand] = kmeans(tempRand', K, 'Replicates', numReplicates, 'Options', opts); 
% %     % Cluster voxels based on 50 singel unit correlation
% %     [IDX_Vox, C, SUMD_Vox] = kmeans(matR_SU, K, 'Replicates', numReplicates, 'Options', opts);
%     
%     
%     Clustering_Rand.resultKMeans(iK).SU_indCluster = IDX_SU_Rand;
%     Clustering_Rand.resultKMeans(iK).SU_sumD = SUMD_SU_Rand;
% %     Clustering.resultKMeans(iK).Vox_indCluster = IDX_Vox;
% %     Clustering.resultKMeans(iK).Vox_sumD = SUMD_Vox;
%     
% end



