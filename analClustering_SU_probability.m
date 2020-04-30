% analyClustering_SU_probability.m
%
% 2017/02/04 SHP
% Compute the probability that one voxel is clustered with another voxel
% 1. Load the "*probability.mat" data 
% 2. For cell-by-cell, mark other voxels that were assigned to the same cluster

clear all;

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjBOLD = 'Art';
dirDataHome = '/procdata/parksh/';
dirDataBOLD= fullfile(dirDataHome, nameSubjBOLD); %fullfile('/procdata/parksh/', nameSubjBOLD); %

% load voxel clustering results
load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  
% load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_voxel_probability_critCorr2.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  
 

for iK=1:length(paramClustering_global.setK)

    
    targetK = paramClustering_global.setK(iK);
    
    fprintf(1, 'K = %d: movie mask \n', targetK);
%     tic;
    matProb = NaN(size(Clustering_moviemask_valid.resultKMeans(1).SU_indCluster, 1));
    for iSU = 1:size(Clustering_moviemask_valid.resultKMeans(1).SU_indCluster, 1)
        matClusterID_SU = repmat(Clustering_moviemask_valid.resultKMeans(iK).SU_indCluster(iSU,:),...
            size(Clustering_moviemask_valid.resultKMeans(iK).SU_indCluster, 1), 1);
        tempMat = [];
        tempMat = Clustering_moviemask_valid.resultKMeans(iK).SU_indCluster - matClusterID_SU;
        tempMat(tempMat~=0) = NaN;
        tempMat(~isnan(tempMat)) = 1;
        tempMat(isnan(tempMat))=0;
        matProb(:,iSU) = sum(tempMat, 2)./100;
    end
    
    Clustering_moviemask_valid.resultKMeans(iK).matProb = matProb;
%     toc;

%     matProb = single(matProb);
%     save(fullfile(dirDataBOLD, sprintf('ClusteringProbability_TorRhoSigSpiArtMovie123_moviemask_%dMeans.mat',targetK)),...
%         'Clustering_moviemask', 'matProb', '-v7.3')
% 
%     fprintf(1, 'K = %d: movie mask: Results saved \n', targetK);

end

save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ...
    'Clustering*', 'param*')  



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

