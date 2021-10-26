% doClusteringVoxTS.m
%

clear all;

% load the data
% nameSubjNeural = 'Spi'; %'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
%
addpath('/library/matlab_utils/')

% Load data files
dirDataHome = '/procdata/parksh/';
% dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% fileNameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

% fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
% load(fullfile(dirDataNeural, filenameNeural))
% load(fullfile(dirDataNeural, fileNameNeural_BLP))
fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
load(fullfile(dirDataBOLD, filenameBOLD))

% fMRI movie-driven activity & brain-only mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');



setMovie = [1 2 3];

% fMRI tc
fmritc=[];
indMovieBOLD = find(ismember(paramBOLD.unimov, setMovie)>0);
for iM = indMovieBOLD %1:length(indMovieBOLD)
    curvoltc = voltcIndMov{iM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc = cat(4,fmritc,pcvoltc(:,:,:,8:125));
end
% Reshape BOLD 4-d data
[nx, ny, nz, nt] = size(fmritc);
nVox = nx*ny*nz;
matVoxTC = reshape(fmritc, nVox, nt);

clear curvoltc avgvoltc pcvoltc fmritc voltcIndMov

%% K-means Clustering based on timeseries

setK = 2:20; %15;
opts = statset('Display','final');
numReplicates = 20; % 100; %5;

paramClustering.methods = 'KMeans';
paramClustering.setK = setK;
paramClustering.numReplicates = numReplicates;
paramClustering.descriptions = 'Clustering of voxels based on the voxel timeseries';


% Apply movie-driven mask to voxel TS matrix
% 2017/01/19 & valid channels
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D
matVoxTC_moviemask = matVoxTC(moviemask_vec,:); % 15495 voxels

brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D
matVoxTC_brainmask = matVoxTC(brainmask_vec,:);  % 27113 voxels

% % totalSS
[a, c, totalSS] = kmeans(matVoxTC_moviemask, 1);
Clustering_moviemask.totalSS = totalSS;
[a, c, totalSS] = kmeans(matVoxTC_brainmask, 1);
Clustering_brainmask.totalSS = totalSS;

% matCluster_SU = NaN(size(matR_SU, 2),
for iK = 1:length(setK)
    
    K = setK(iK);
    
    fprintf(1, ':: K = %d; Movie-driven mask ::\n', K);
    % Cluster voxels based on 50 singel unit correlation
    [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matVoxTC_moviemask, K, 'Replicates', numReplicates, 'Options', opts);
    
    Clustering_moviemask.resultKMeans(iK).Vox_indCluster = IDX_Vox;
    Clustering_moviemask.resultKMeans(iK).Vox_sumD = SUMD_Vox;
    
    fprintf(1, ':: K = %d; brain mask ::\n', K);
    % Cluster voxels based on 50 singel unit correlation
    [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matVoxTC_brainmask, K, 'Replicates', numReplicates, 'Options', opts);
    
    Clustering_brainmask.resultKMeans(iK).Vox_indCluster = IDX_Vox;
    Clustering_brainmask.resultKMeans(iK).Vox_sumD = SUMD_Vox;
    
end


save(fullfile(dirDataBOLD, sprintf('ClusteringVoxTS_%sMovie123.mat', nameSubjBOLD)), 'Clustering*', 'paramClustering')

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iK=1:length(setK)
    plot(iK+1, mean(Clustering_moviemask.resultKMeans(iK).Vox_sumD), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
end
xlim([1 20])
set(gca, 'LineWidth', 2, 'FontSize', 12)
box off
title('Clustering of voxels')
xlabel('Number of cluster')
ylabel('Within-cluster distance')

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
%     [IDX_SU, C, SUMD_SU] = kmeans(matFR_SU', K, 'Replicates', numReplicates, 'Options', opts);
%     % Cluster time points based on 50 singel unit correlation
%     [IDX_t, C, SUMD_t] = kmeans(matFR_SU, K, 'Replicates', numReplicates, 'Options', opts);
%
%     % Cluster single units based on normalized (z-scored) time series
%     [IDX_SUnorm, C, SUMD_SUnorm] = kmeans(matFR_SU_norm', K, 'Replicates', numReplicates, 'Options', opts);
%
%
%     ClusteringSDF.resultKMeans(iK).SU_indCluster = IDX_SU;
%     ClusteringSDF.resultKMeans(iK).SU_sumD = SUMD_SU;
%     ClusteringSDF.resultKMeans(iK).time_indCluster = IDX_t;
%     ClusteringSDF.resultKMeans(iK).time_sumD = SUMD_t;
%
%     ClusteringSDFnorm.resultKMeans(iK).SU_indCluster = IDX_SUnorm;
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

