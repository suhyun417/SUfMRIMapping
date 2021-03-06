% analyClustering_SU_probability.m
%
% 2021/06/21 SHP
% Apply this to multiple face patch clustering results using only face cells
% 2020/05/03 SHP
% Apply this to multiple face patch clustering results
% 2017/02/04 SHP
% Compute the probability that one cell is clustered with another cell
% 1. Load the "*probability.mat" data 
% 2. For cell-by-cell, mark other voxels that were assigned to the same cluster

clear all;

% setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjBOLD = 'Art';
dirDataHome = '/procdata/parksh/_macaque';
dirDataBOLD= fullfile(dirDataHome, nameSubjBOLD); %fullfile('/procdata/parksh/', nameSubjBOLD); %

% load  clustering results
% fname = fullfile(dirDataBOLD, 'Clustering_CorrMap_4FPs_Movie123_probability.mat'); %, 'Clustering_brainmask', 'param*') 
% fname = fullfile(dirDataBOLD, 'Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat');

fname = fullfile(dirDataBOLD, 'Clustering_CorrMap_4FPs_faceselective_Movie123_probability.mat'); %, 'Clustering_brainmask', 'param*') 
load(fname);
% load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  
% load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_voxel_probability_critCorr2.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  
 
Clustering = Clustering_brainmask; %Clustering_meanROI; %Clustering_brainmask; %Clustering_meanROI; %Clustering_moviemask_valid;
paramClustering_global.numRepeat = 100;

for iK=1:length(paramClustering_global.setK)

    
    targetK = paramClustering_global.setK(iK);
    
%     fprintf(1, 'K = %d: movie mask \n', targetK);
%     tic;
    matProb = NaN(size(Clustering.resultKMeans(1).SU_indCluster, 1));
    for iSU = 1:size(Clustering.resultKMeans(1).SU_indCluster, 1)
        matClusterID_SU = repmat(Clustering.resultKMeans(iK).SU_indCluster(iSU,:),...
            size(Clustering.resultKMeans(iK).SU_indCluster, 1), 1);
        tempMat = [];
        tempMat = Clustering.resultKMeans(iK).SU_indCluster - matClusterID_SU; %
        tempMat(tempMat~=0) = NaN; % cells that are not in the same cluster changed to NaN
        tempMat(~isnan(tempMat)) = 1; % same cluster changed to 1
        tempMat(isnan(tempMat))=0; % different cluster changed to 0
        matProb(:,iSU) = sum(tempMat, 2)./paramClustering_global.numRepeat; % sum across repeats
    end
    
    Clustering.resultKMeans(iK).matProb = matProb;
%     toc;

%     matProb = single(matProb);
%     save(fullfile(dirDataBOLD, sprintf('ClusteringProbability_TorRhoSigSpiArtMovie123_moviemask_%dMeans.mat',targetK)),...
%         'Clustering_moviemask', 'matProb', '-v7.3')
% 
%     fprintf(1, 'K = %d: movie mask: Results saved \n', targetK);

end

Clustering_brainmask = Clustering;
% Clustering_meanROI = Clustering; %
clear Clustering

save(fname, 'Clustering*', 'param*', '-append');

% save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ...
%     'Clustering*', 'param*')  



% TEMP
matProb = NaN(size(SU_indCluster, 1));
    for iSU = 1:size(SU_indCluster, 1)
        matClusterID_SU = repmat(SU_indCluster(iSU,:),...
            size(SU_indCluster, 1), 1);
        tempMat = [];
        tempMat = SU_indCluster - matClusterID_SU;
        tempMat(tempMat~=0) = NaN;
        tempMat(~isnan(tempMat)) = 1;
        tempMat(isnan(tempMat))=0;
        matProb(:,iSU) = sum(tempMat, 2)./100;
    end

