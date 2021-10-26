function [] = doClusteringCorrMap()
% 2017/01/23 SHP modified from doClusteringCorrMap.m
% Apply K-means clustering to correlation map of single units to categorize
% cells based on their correlation pattern with other brain areas
% Now apply confidence interval computed from separate bootstrap procedure,
% use voxels that were significant (i.e. correlation value is outside of
% the 95% confidence interval)
% 2016/06/19 Add clustering on masked map (movie-driven mask applied)

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
nameSubjNeural = 'Tor'; %'Spi'; %'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

%% Load the corr map
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');

% load the masks (movie-driven & brain-only mask)
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');


% %% K-means Clustering based on correlation maps
% 
% setK = 2:15;
% opts = statset('Display','final');
% numReplicates = 30; %5;
% 
% % 2017/01/19: Excluded some units, in case of Spice's data which hasn't been went through this procedure yet 
% switch lower(nameSubjNeural)
%     case 'spi'
%         excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
%         validChanIndex_clustering = setdiff(paramCorr.validChanIndex, excChanIndex);
%         validChanID_clustering = paramCorr.validChanID(validChanIndex_clustering,:);
%     otherwise
%         validChanIndex_clustering = 1:length(paramCorr.validChanIndex);
%         validChanID_clustering = paramCorr.validChanID;
% end 
% 
% paramClustering.validChanIndex = validChanIndex_clustering;
% paramClustering.validChanID = validChanID_clustering;
% 
% Clustering.methods = 'KMeans';
% Clustering.setK = setK;
% Clustering.numReplicates = numReplicates;
% Clustering.descriptions = 'Clustering of neurons based on the whole-brain correlation maps';
% 
% 
% % totalSS
% [a, c, totalSS] = kmeans(matR_SU(:, validChanIndex_clustering)', 1);
% Clustering.totalSS = totalSS;
% 
% diary(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_ReplicatesRecords', nameSubjNeural, nameSubjBOLD)))
% for iK = 1:length(setK)
%     
%     K = setK(iK);
%     
%     diary on
%     fprintf(1, ':: K = %d; ::\n', K);
%     % Cluster single units based on whole brain correlation
%     [IDX_SU, C_SU, SUMD_SU, totalD] = kmeans(matR_SU', K, 'Replicates', numReplicates, 'Options', opts); 
%     
%     diary off
%     % Cluster voxels based on 50 singel unit correlation
%     [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU, K, 'Replicates', numReplicates, 'Options', opts);
% 
%     
%     Clustering.resultKMeans(iK).SU_indCluster = IDX_SU;
%     Clustering.resultKMeans(iK).SU_sumD = SUMD_SU;
%     Clustering.resultKMeans(iK).SU_totalD = totalD;
%     Clustering.resultKMeans(iK).Vox_indCluster = IDX_Vox;
%     Clustering.resultKMeans(iK).Vox_sumD = SUMD_Vox;
%     
% end
% save(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'Clustering');
% 
% %%%   
% % plot Clustering results
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% for iK=1:length(setK)
%     plot(iK+1, Clustering.resultKMeans(iK).SU_sumD, 'ko', 'MarkerSize', 10, 'LineWidth', 2)
%     hold on
% end
% xlim([1 16])
% set(gca, 'LineWidth', 2, 'FontSize', 12)
% box off
% title('Clustering of single units')
% xlabel('Number of cluster')
% ylabel('Within-cluster distance')
% % print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_SU'), '-depsc');
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
% 

%% 2. Clustering on masked map 

setK = 2:15;
opts = statset('Display','final');
numReplicates = 100; % 5; %10; %5;

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
paramClustering.descriptions = 'Clustering of neurons based on the masked whole-brain correlation maps';

%% Apply confidence intervals computed from separate bootstrap procedure
dirDataCI = fullfile(dirDataNeural, 'corrMap_resultsCI/resultsCI_eachCell');
sigVox_CI95=[];
for iUnit = 1:length(validChanIndex_clustering)    
    
    cellID = validChanID_clustering(iUnit,:);
    fprintf(1, 'Unit # %d, Cell ID: %s \n', iUnit, cellID);
   
    fileNameCI = sprintf('corrcoeffCI_%s%s_cell%d.mat',...
        nameSubjNeural, nameSubjBOLD, iUnit);
    load(fullfile(dirDataCI, fileNameCI))
    
    catCI95 = cat(1, resultBS(1).VoxCI.matCI95);
%     catCI99 = cat(1, resultBS(1).VoxCI.matCI99);
    
%     CI95_in = catCI95(mask_amp1_sub, :);
%     setCI95_in(iUnit, :) = mean(CI95_in);
%     % catRho_org = cat(1, resultBS(1).VoxCI.rho_org).*(-1);
    corrMap_org = matR_SU(:,iUnit);

    aa=cat(2, corrMap_org<catCI95(:,1), corrMap_org>catCI95(:,2));
    sigVox_CI95(:,iUnit) = sum(aa,2);    
    
    matR_SU(sum(aa,2)<1, iUnit) = 0; % replace non-significant
%     correlations to zero (mean of all those non-significant voxels were
%     0.0075, median is 1.18e-04)
end
    
 



%%
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

% 2017/01/19 Record the results of each replicate (number of iterations,
% sum of the distance)
% diary(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked_ReplicatesRecords', nameSubjNeural, nameSubjBOLD)))
for iK = 1:length(setK)
    
    K = setK(iK);                
    
%     diary on
    fprintf(1, ':: K = %d; Movie-driven mask ::\n', K);
    % 1. Movie-driven mask
    % Cluster single units based on whole brain correlation
    [IDX_SU, C_SU, SUMD_SU, totalD] = kmeans(matR_SU_moviemask', K, 'Replicates', numReplicates, 'Options', opts); 
    
%     diary off 
    % Cluster voxels based on 50 singel unit correlation
    [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_moviemask, K, 'Replicates', 5, 'Options', opts); % numReplicates, 'Options', opts);
    
    Clustering_moviemask.resultKMeans(iK).SU_indCluster = IDX_SU;
    Clustering_moviemask.resultKMeans(iK).SU_sumD = SUMD_SU;
    Clustering_moviemask.resultKMeans(iK).SU_totalD = totalD;
    Clustering_moviemask.resultKMeans(iK).SU_centeroid = C_SU;
    Clustering_moviemask.resultKMeans(iK).Vox_indCluster = IDX_Vox;
    Clustering_moviemask.resultKMeans(iK).Vox_sumD = SUMD_Vox;

    
%     diary on
    fprintf(1, ':: K = %d; Brain mask ::\n', K);
    % 2. Brain mask
    % Cluster single units based on whole brain correlation
    [IDX_SU, C_SU, SUMD_SU, totalD] = kmeans(matR_SU_brainmask', K, 'Replicates', numReplicates, 'Options', opts); 
    
%     diary off
    % Cluster voxels based on 50 singel unit correlation
    [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_brainmask, K, 'Replicates', 5, 'Options', opts); %, numReplicates, 'Options', opts);
    
    Clustering_brainmask.resultKMeans(iK).SU_indCluster = IDX_SU;
    Clustering_brainmask.resultKMeans(iK).SU_sumD = SUMD_SU;
    Clustering_brainmask.resultKMeans(iK).SU_totalD = totalD;
    Clustering_brainmask.resultKMeans(iK).SU_centeroid = C_SU;
    Clustering_brainmask.resultKMeans(iK).Vox_indCluster = IDX_Vox;
    Clustering_brainmask.resultKMeans(iK).Vox_sumD = SUMD_Vox;
    
end
save(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked_significant.mat', nameSubjNeural, nameSubjBOLD)),...
    'Clustering_moviemask', 'Clustering_brainmask', 'paramClustering');

%%%   
% plot Clustering results
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iK=1:length(setK)
    plot(iK+1, mean(Clustering_moviemask.resultKMeans(iK).SU_sumD), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
end
xlim([1 16])
set(gca, 'LineWidth', 2, 'FontSize', 12)
box off
title('Clustering of single units')
xlabel('Number of cluster')
ylabel('Within-cluster distance')
% print(gcf, fullfile(dirFig, 'kmeans_distanceElbowPlot_SU'), '-depsc');


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



