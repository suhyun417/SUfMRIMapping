function [] = doClusteringCorrMap()


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
nameSubjNeural = 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

%% Load the corr map
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');

% 4) Movie-driven mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');

[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;

% Apply movie-driven mask to correlation matrix
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D
matR_SU_moviemask = matR_SU(moviemask_vec,:); % 15495 voxels

%% Gaussian Mixture model
setK = 1:15;

AIC = zeros(1,length(setK));
GMModels = cell(1,length(setK));
options = statset('MaxIter',500,'Display','final');

for k = setK %1:7 %setK %1:4
    GMModels{k} = fitgmdist(matR_SU_moviemask,k,'Options',options,'CovType','diagonal'); %fitgmdist(matR_SU,k,'Options',options,'CovType','diagonal');
    AIC(k)= GMModels{k}.AIC;
end

[minAIC,numComponents] = min(AIC);
numComponents

BestModel = GMModels{numComponents}


% 
% 
% %% K-means Clustering based on correlation maps
% 
% setK = 2:15;
% opts = statset('Display','final');
% numReplicates = 5;
% 
% Clustering.methods = 'KMeans';
% Clustering.setK = setK;
% Clustering.numReplicates = numReplicates;
% 
% % totalSS
% [a, c, totalSS] = kmeans(matR_SU', 1);
% Clustering.totalSS = totalSS;
% 
% % matCluster_SU = NaN(size(matR_SU, 2), 
% for iK = 1:length(setK)
%     
%     K = setK(iK);
%     % Cluster single units based on whole brain correlation
%     [IDX_SU, C, SUMD_SU, totalD] = kmeans(matR_SU', K, 'Replicates', numReplicates, 'Options', opts); 
%     % Cluster voxels based on 50 singel unit correlation
%     [IDX_Vox, C, SUMD_Vox] = kmeans(matR_SU, K, 'Replicates', numReplicates, 'Options', opts);
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



