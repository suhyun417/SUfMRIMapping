function [] = doClusteringCorrMap_voxel_probability()
% 2017/01/21 SHP: modified from doClusteringCorrMap.m
% Repeat K-means clustering multiple times to categorize
% cells based on their correlation pattern with other brain areas,
% save the results of K-means each time, then compute the probability of
% each cell is clusterd as same cluster with other cells
% Use only the maps masked with movie-driven masks

clear all;

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
nameSubjNeural = 'Tor'; %'Spi'; %'Tor'; %'Spi'; %'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

%% Load the corr map
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');

% load the masks (movie-driven & brain-only mask)
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');


%% 2. Clustering on masked map 
numRepeat = 50; %100; % number of repetition for entire clustering

setK = 2:20; %2:15;
opts = statset('Display','final');
numReplicates = 5; %30; %100; % 5; %10; %5;

% Excluded some units 
switch lower(nameSubjNeural)
    case 'spi'
        excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
        validChanIndex_clustering = setdiff(paramCorr.validChanIndex, excChanIndex);
        validChanID_clustering = paramCorr.validChanID(validChanIndex_clustering,:);
    otherwise
        validChanIndex_clustering = 1:length(paramCorr.validChanIndex);
        validChanID_clustering = paramCorr.validChanID;
end 

paramClustering.validChanIndex = validChanIndex_clustering;
paramClustering.validChanID = validChanID_clustering;

paramClustering.methods = 'KMeans';
paramClustering.setK = setK;
paramClustering.numReplicates = numReplicates;
paramClustering.descriptions = 'Clustering of neurons based on the masked whole-brain correlation maps: multiple repetitions (numRepeat) of an execution of kmeans function, which had its own multiple "replicates" (numReplicates)';
paramClustering.numRepeat = numRepeat;



[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;

% Apply movie-driven mask to correlation matrix 
% 2017/01/19 & valid channels 
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D
matR_SU_moviemask = matR_SU(moviemask_vec,validChanIndex_clustering); %matR_SU(moviemask_vec,:); % 15495 voxels

brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
matR_SU_brainmask = matR_SU(brainmask_vec,validChanIndex_clustering); %matR_SU(brainmask_vec,:); % 27113 voxels

% % totalSS
[a, c, totalSS] = kmeans(matR_SU_moviemask', 1);
Clustering_moviemask.totalSS = totalSS;
[a, c, totalSS] = kmeans(matR_SU_brainmask', 1);
Clustering_brainmask.totalSS = totalSS;

% % 2017/01/19 Record the results of each replicate (number of iterations,
% % sum of the distance)
% diary(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked_ReplicatesRecords', nameSubjNeural, nameSubjBOLD)))
for iK = 1:length(setK)
    
    K = setK(iK);                
    
    Vox_indCluster_moviemask = NaN(sum(moviemask_vec), numRepeat);
    Vox_sumD_moviemask = NaN(K, numRepeat);
    
    Vox_indCluster_brainmask = NaN(sum(brainmask_vec), numRepeat);
    Vox_sumD_brainmask = NaN(K, numRepeat);

    for iRep = 1:numRepeat
        
        %     diary on
        fprintf(1, ':: K = %d; Movie-driven mask ::\n', K);
        % 1. Movie-driven mask
%         % Cluster single units based on whole brain correlation
%         [IDX_SU, C_SU, SUMD_SU, totalD] = kmeans(matR_SU_moviemask', K, 'Replicates', numReplicates, 'Options', opts);
%         
%         SU_indCluster_moviemask(:, iRep) = IDX_SU;
%         SU_sumD_moviemask(:, iRep) = SUMD_SU;
        
        %     diary off
            % Cluster voxels based on 50 singel unit correlation
            [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_moviemask, K, 'Replicates', numReplicates, 'Options', opts);
            
            Vox_indCluster_moviemask(:, iRep) = IDX_Vox;
            Vox_sumD_moviemask(:, iRep) = SUMD_Vox;
        
        %     Clustering_moviemask.resultKMeans(iK).SU_indCluster = IDX_SU;
        %     Clustering_moviemask.resultKMeans(iK).SU_sumD = SUMD_SU;
        %     Clustering_moviemask.resultKMeans(iK).SU_totalD = totalD;
        %     Clustering_moviemask.resultKMeans(iK).SU_centeroid = C_SU;
%             Clustering_moviemask.resultKMeans(iK).Vox_indCluster = IDX_Vox;
%             Clustering_moviemask.resultKMeans(iK).Vox_sumD = SUMD_Vox;
        
        %     diary on
        fprintf(1, ':: K = %d; Brain mask ::\n', K);
        % 2. Brain mask
%         % Cluster single units based on whole brain correlation
%         [IDX_SU, C_SU, SUMD_SU, totalD] = kmeans(matR_SU_brainmask', K, 'Replicates', numReplicates, 'Options', opts);
%         
%         SU_indCluster_brainmask(:, iRep) = IDX_SU;
%         SU_sumD_brainmask(:, iRep) = SUMD_SU;
        
        %     diary off
            % Cluster voxels based on 50 singel unit correlation
            [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_brainmask, K, 'Replicates', numReplicates, 'Options', opts);
        
            Vox_indCluster_brainmask(:, iRep) = IDX_Vox;
            Vox_sumD_brainmask(:, iRep) = SUMD_Vox;
            
        %     Clustering_brainmask.resultKMeans(iK).SU_indCluster = IDX_SU;
        %     Clustering_brainmask.resultKMeans(iK).SU_sumD = SUMD_SU;
        %     Clustering_brainmask.resultKMeans(iK).SU_totalD = totalD;
        %     Clustering_brainmask.resultKMeans(iK).SU_centeroid = C_SU;
        %     Clustering_brainmask.resultKMeans(iK).Vox_indCluster = IDX_Vox;
        %     Clustering_brainmask.resultKMeans(iK).Vox_sumD = SUMD_Vox;
    end
    
    Clustering_moviemask.resultKMeans(iK).Vox_indCluster = Vox_indCluster_moviemask;
    Clustering_moviemask.resultKMeans(iK).Vox_sumD = Vox_sumD_moviemask;
    
    Clustering_brainmask.resultKMeans(iK).Vox_indCluster = Vox_indCluster_brainmask;
    Clustering_brainmask.resultKMeans(iK).Vox_sumD = Vox_sumD_brainmask;
    
end
save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_voxels_1.mat', nameSubjNeural, nameSubjBOLD)),...
    'Clustering_moviemask', 'Clustering_brainmask', 'paramClustering');

%%%   
% plot Clustering results
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iK=1:length(setK)
    plot(iK+1, mean(Clustering_moviemask.resultKMeans(iK).Vox_sumD), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
end
xlim([1 16])
set(gca, 'LineWidth', 2, 'FontSize', 12)
box off
title('Clustering of voxels')
xlabel('Number of cluster')
ylabel('Within-cluster distance')
% print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_SU'), '-depsc');