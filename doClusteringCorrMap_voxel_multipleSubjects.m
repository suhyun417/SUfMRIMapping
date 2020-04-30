function [] = doClusteringCorrMap_voxel_multipleSubjects()
% 2017/01/22 SHP
% Apply K-means clustering to correlation map of single units from multiple
% subjects
% Modified from 'doClusteringCorrMap.m'

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
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% load the masks (movie-driven & brain-only mask)
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');


%% Load the corr map and select valid channel: subject by subject
numSubject = size(setNameSubjNeural, 2);
matR_SU_all = [];

for iSubj = 1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');
    
    % 2017/01/19: Excluded some units, in case of Spice's data which hasn't been went through this procedure yet
    switch lower(nameSubjNeural)
        case 'spi'
            excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
            validChanIndex_clustering = setdiff(paramCorr.validChanIndex, excChanIndex);
            validChanID_clustering = cat(2, paramCorr.validChanID(validChanIndex_clustering,:), num2str(ones(size(validChanIndex_clustering)).*iSubj)); %paramCorr.validChanID(validChanIndex_clustering,:);
        otherwise
            validChanIndex_clustering = (1:length(paramCorr.validChanIndex))';
            validChanID_clustering = cat(2, paramCorr.validChanID, num2str(ones(size(paramCorr.validChanIndex)).*iSubj)); %paramCorr.validChanID;
    end
    
    paramClustering(iSubj).nameSubj = nameSubjNeural;
    paramClustering(iSubj).validChanIndex = validChanIndex_clustering;
    paramClustering(iSubj).validChanID = validChanID_clustering;
    
    matR_SU_valid = matR_SU(:, validChanIndex_clustering);
    matR_SU_all = cat(2, matR_SU_all, matR_SU_valid);
    
     clear matR_SU matR_SU_valid
end

% %% K-means Clustering based on correlation maps
% 
% setK = 2:20; %15;
% opts = statset('Display','final');
% numReplicates = 10; %30; %5;
% 
% Clustering.methods = 'KMeans';
% Clustering.setK = setK;
% Clustering.numReplicates = numReplicates;
% Clustering.descriptions = 'Clustering of neurons based on the whole-brain correlation maps';
% 
% 
% % totalSS
% [a, c, totalSS] = kmeans(matR_SU_all, 1);
% Clustering.totalSS = totalSS;
% 
% % diary(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_ReplicatesRecords', nameSubjNeural, nameSubjBOLD)))
% for iK = 1:length(setK)
%     
%     K = setK(iK);
%     
% %     diary on
%     fprintf(1, ':: K = %d; ::\n', K);
% %     % Cluster single units based on whole brain correlation
% %     [IDX_SU, C_SU, SUMD_SU, totalD] = kmeans(matR_SU_all', K, 'Replicates', numReplicates, 'Options', opts); 
%     
% %     diary off
%     % Cluster voxels based on 50 singel unit correlation
%     [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_all, K, 'Replicates', numReplicates, 'Options', opts);
% 
%     
% %     Clustering.resultKMeans(iK).SU_indCluster = IDX_SU;
% %     Clustering.resultKMeans(iK).SU_sumD = SUMD_SU;
% %     Clustering.resultKMeans(iK).SU_totalD = totalD;
%     Clustering.resultKMeans(iK).Vox_indCluster = IDX_Vox;
%     Clustering.resultKMeans(iK).Vox_sumD = SUMD_Vox;
%     
% end
% save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_voxel.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), 'Clustering'); %nameSubjNeural, nameSubjBOLD)), 'Clustering');
% 
% %%%   
% % plot Clustering results
% % figure;
% % set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% % for iK=1:length(setK)
% %     plot(iK+1, Clustering.resultKMeans(iK).SU_sumD, 'ko', 'MarkerSize', 10, 'LineWidth', 2)
% %     hold on
% % end
% % xlim([1 16])
% % set(gca, 'LineWidth', 2, 'FontSize', 12)
% % box off
% % title('Clustering of single units')
% % xlabel('Number of cluster')
% % ylabel('Within-cluster distance')
% % % print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_SU'), '-depsc');
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% for iK=1:length(setK)
%     plot(iK+1, Clustering.resultKMeans(iK).Vox_sumD, 'ko', 'MarkerSize', 10, 'LineWidth', 2)
%     hold on
% end
% xlim([1 16])
% set(gca, 'LineWidth', 2, 'FontSize', 12)
% box off
% title('Clustering of voxels')
% xlabel('Number of cluster')
% ylabel('Within-cluster distance')
% % print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_Vox'), '-depsc');


%% 2. Clustering on masked map 

setK = 2:20; %15;
opts = statset('Display','final');
numReplicates = 100; %10; %20; %30; %100; % 5; %10; %5;


paramClustering_global.methods = 'KMeans';
paramClustering_global.setK = setK;
paramClustering_global.numReplicates = numReplicates;
paramClustering_global.descriptions = 'Clustering of neurons based on the masked whole-brain correlation maps';

if flagParallel
    pool = parpool;                      % Invokes workers
    stream = RandStream('mlfg6331_64');  % Random number stream
    opts = statset('UseParallel',1,'UseSubstreams',1, 'Streams',stream,...
        'MaxIter', 1000, 'Display','final');
    paramClustering_global.parallel = opts;
end

[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;

% Apply movie-driven mask to correlation matrix 
% 2017/01/19 & valid channels 
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D
matR_SU_all_moviemask = matR_SU_all(moviemask_vec,:); %matR_SU(moviemask_vec,:); % 15495 voxels

brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
matR_SU_all_brainmask = matR_SU_all(brainmask_vec,:); %matR_SU(brainmask_vec,:); % 27113 voxels

% % totalSS
[a, c, totalSS] = kmeans(matR_SU_all_moviemask, 1);
Clustering_moviemask.totalSS = totalSS;
[a, c, totalSS] = kmeans(matR_SU_all_brainmask, 1);
Clustering_brainmask.totalSS = totalSS;

% 2017/01/19 Record the results of each replicate (number of iterations,
% sum of the distance)
% diary(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked_ReplicatesRecords', nameSubjNeural, nameSubjBOLD)))
for iK = 1:length(setK)
    
    K = setK(iK);                
    
%     diary on
    fprintf(1, ':: K = %d; Movie-driven mask ::\n', K);
%     % 1. Movie-driven mask
%     % Cluster single units based on whole brain correlation
%     [IDX_SU, C_SU, SUMD_SU, totalD] = kmeans(matR_SU_all_moviemask', K, 'Replicates', numReplicates, 'Options', opts); 
    
%     diary off 
    % Cluster voxels based on 50 singel unit correlation
    [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_all_moviemask, K, 'Replicates', numReplicates, 'Options', opts);
    
%     Clustering_moviemask.resultKMeans(iK).SU_indCluster = IDX_SU;
%     Clustering_moviemask.resultKMeans(iK).SU_sumD = SUMD_SU;
%     Clustering_moviemask.resultKMeans(iK).SU_totalD = totalD;
%     Clustering_moviemask.resultKMeans(iK).SU_centeroid = C_SU;
    Clustering_moviemask.resultKMeans(iK).Vox_indCluster = IDX_Vox;
    Clustering_moviemask.resultKMeans(iK).Vox_sumD = SUMD_Vox;

end

for iK = 1:length(setK)
    
    K = setK(iK);             
    
%     diary on
    fprintf(1, ':: K = %d; Brain mask ::\n', K);
    % 2. Brain mask
%     % Cluster single units based on whole brain correlation
%     [IDX_SU, C_SU, SUMD_SU, totalD] = kmeans(matR_SU_all_brainmask', K, 'Replicates', numReplicates, 'Options', opts); 
    
%     diary off
    % Cluster voxels based on 50 singel unit correlation
    [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_all_brainmask, K, 'Replicates', 5, 'Options', opts); %, numReplicates, 'Options', opts);
    
%     Clustering_brainmask.resultKMeans(iK).SU_indCluster = IDX_SU;
%     Clustering_brainmask.resultKMeans(iK).SU_sumD = SUMD_SU;
%     Clustering_brainmask.resultKMeans(iK).SU_totalD = totalD;
%     Clustering_brainmask.resultKMeans(iK).SU_centeroid = C_SU;
    Clustering_brainmask.resultKMeans(iK).Vox_indCluster = IDX_Vox;
    Clustering_brainmask.resultKMeans(iK).Vox_sumD = SUMD_Vox;
    
end
save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_voxel.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
    'Clustering_moviemask', 'Clustering_brainmask', 'paramClustering*');

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
% print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_Vox'), '-depsc');


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



