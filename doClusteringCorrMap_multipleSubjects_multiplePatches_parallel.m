function [] = doClusteringCorrMap_multipleSubjects_multiplePatches_parallel(curK, flagSave)
%
% 2020/06/16 SHP 
%       modified to cluster only face-selective neurons from 4
%       face patches. removed the parts using only movie-driven voxels
% 2020/08/13 SHP
%       K-means clustering of correlation maps of single units from 4 face
%       patches using biowulf


% Set directories
dirDataHome =  '/data/parks20/procdata/NeuroMRI'; % Biowulf '/procdata/parksh/'; %nifstorage %'/data/parks20/procdata/NeuroMRI'; % Biowulf
% dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, 'Art');
dirSaveFile = fullfile(dirDataBOLD, 'clustering_multipleFP'); % fullfile(dirDataNeural, 'corrMap_resultsCI_10hz');
if ~exist(dirSaveFile, 'dir')
    mkdir(dirSaveFile)
end


%% Load the data
fname = 'matR4clusteringmultiplepatches_faceselective.mat'; %'matR4clusteringmultiplepatches.mat';
load(fullfile(dirDataBOLD, fname))

% Take care of possible string-double issue from compile
if ischar(curK); eval(sprintf('curK = %s;', curK)); end


%% Clustering on masked map
numRepeat = 100; % number of repetition for entire clustering

setK = curK; %stK:edK;
opts = statset('Display','final');
numReplicates = 10; %5; %100;

paramClustering_global.methods = 'KMeans';
paramClustering_global.setK = setK;
paramClustering_global.numReplicates = numReplicates;
paramClustering_global.descriptions = 'Clustering of neurons based on the whole-brain correlation maps with brainmask and moviemask';
paramClustering_global.numRepeat = numRepeat;



% totalSS: brainmask
[a, c, totalSS] = kmeans(matR_brainmask', 1); % single unit
Clustering_brainmask.totalSS_SU = totalSS;
[a, c, totalSS] = kmeans(matR_brainmask, 1); % voxel
Clustering_brainmask.totalSS_Vox = totalSS;

% % totalSS: moviemask
% [a, c, totalSS] = kmeans(matR_moviemask', 1); % single unit
% Clustering_moviemask.totalSS_SU = totalSS;
% [a, c, totalSS] = kmeans(matR_moviemask, 1); % voxel
% Clustering_moviemask.totalSS_Vox = totalSS;

% for K = setK
    
%     K = setK(iK);
        
    %% brainmask
    SU_indCluster_brainmask = NaN(size(matR_brainmask, 2), numRepeat);
    SU_sumD_brainmask = NaN(curK, numRepeat);
    
    Vox_indCluster_brainmask = NaN(size(matR_brainmask, 1), numRepeat);
    Vox_sumD_brainmask = NaN(curK, numRepeat);
    
    for iRep = 1:numRepeat
        
        fprintf(1, ':: K = %d; Rep = %d/%d; brainmask ::\n', curK, iRep, numRepeat);
        
        % Cluster single units based on whole brain correlation
        [IDX_SU, C_SU, SUMD_SU] = kmeans(matR_brainmask', curK,...
            'Replicates', numReplicates, 'Options', opts);
        
        SU_indCluster_brainmask(:, iRep) = IDX_SU;
        SU_sumD_brainmask(:, iRep) = SUMD_SU;
        
        % Cluster voxels based on single unit correlation
        [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_brainmask, curK,...
            'Replicates', numReplicates, 'Options', opts);
        
        Vox_indCluster_brainmask(:, iRep) = IDX_Vox;
        Vox_sumD_brainmask(:, iRep) = SUMD_Vox;
        
    end
    
    Clustering_brainmask.resultKMeans.SU_indCluster = SU_indCluster_brainmask;
    Clustering_brainmask.resultKMeans.SU_sumD = SU_sumD_brainmask;
    Clustering_brainmask.resultKMeans.Vox_indCluster = Vox_indCluster_brainmask;
    Clustering_brainmask.resultKMeans.Vox_sumD = Vox_sumD_brainmask;
    
    if flagSave
        save(fullfile(dirSaveFile, sprintf('Clustering_CorrMap_4FPs_faceselective_Movie123_probability_K%d.mat', curK)),...
            'Clustering*', 'paramClustering*');
%         save(fullfile(dirSaveFile, sprintf('Clustering_CorrMap_4FPs_Movie123_probability_K%d.mat', curK)),...
%             'Clustering*', 'paramClustering*');
% %         save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
% %             'Clustering_moviemask', 'paramClustering*');
        fprintf(1, ':: K = %d; Results saved \n', curK);
    end
    
%     %% moviemask
%     SU_indCluster_moviemask = NaN(size(matR_moviemask, 2), numRepeat);
%     SU_sumD_moviemask = NaN(curK, numRepeat);
%     
%     Vox_indCluster_moviemask = NaN(size(matR_moviemask, 1), numRepeat);
%     Vox_sumD_moviemask = NaN(curK, numRepeat);
%     
%     for iRep = 1:numRepeat
%         
%         fprintf(1, ':: K = %d; Rep = %d/%d; moviemask ::\n', curK, iRep, numRepeat);
%         
%         % Cluster single units based on whole brain correlation
%         [IDX_SU, C_SU, SUMD_SU] = kmeans(matR_moviemask', curK,...
%             'Replicates', numReplicates, 'Options', opts);
%         
%         SU_indCluster_moviemask(:, iRep) = IDX_SU;
%         SU_sumD_moviemask(:, iRep) = SUMD_SU;
%         
%         % Cluster voxels based on single unit correlation
%         [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_moviemask, curK,...
%             'Replicates', numReplicates, 'Options', opts);
%         
%         Vox_indCluster_moviemask(:, iRep) = IDX_Vox;
%         Vox_sumD_moviemask(:, iRep) = SUMD_Vox;
%         
%     end
%     
%     Clustering_moviemask.resultKMeans.SU_indCluster = SU_indCluster_moviemask;
%     Clustering_moviemask.resultKMeans.SU_sumD = SU_sumD_moviemask;
%     Clustering_moviemask.resultKMeans.Vox_indCluster = Vox_indCluster_moviemask;
%     Clustering_moviemask.resultKMeans.Vox_sumD = Vox_sumD_moviemask;
    
%     if flagSave
%         save(fullfile(dirSaveFile, sprintf('Clustering_CorrMap_4FPs_Movie123_probability_K%d.mat', curK)),...
%             'Clustering*', 'paramClustering*');
% %         save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
% %             'Clustering_moviemask', 'paramClustering*');
%         fprintf(1, ':: K = %d; Results saved \n', curK);
%     end
    
% end


%% movie-mask
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



