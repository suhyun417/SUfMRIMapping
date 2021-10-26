% analyClustering_voxel_probability.m
%
% 2017/01/28 SHP
% Compute the probability that one voxel is clustered with another voxel
% 1. Load the "*probability.mat" data 
% 2. For voxel-by-voxel, mark other voxels that were assigned to the same cluster

clear all;

nameSubjBOLD = 'Art';
dirDataBOLD= fullfile('/data/parks20/procdata/NeuroMRI/', nameSubjBOLD); %fullfile('/procdata/parksh/', nameSubjBOLD); %

load(fullfile(dirDataBOLD, 'Clustering_TorRhoSigSpiArtMovie123_new_masked_voxel_probability.mat'))

% try
pool=parpool;
% end

for iK=1:length(paramClustering_global.setK)

    
    targetK = paramClustering_global.setK(iK);
    
    fprintf(1, 'K = %d: movie mask \n', targetK);
%     tic;
    matProb = NaN(size(Clustering_moviemask.resultKMeans(1).Vox_indCluster, 1));
    parfor iVox = 1:size(Clustering_moviemask.resultKMeans(1).Vox_indCluster, 1)
        matClusterID_vox = repmat(Clustering_moviemask.resultKMeans(iK).Vox_indCluster(iVox,:),...
            size(Clustering_moviemask.resultKMeans(iK).Vox_indCluster, 1), 1);
        tempMat = [];
        tempMat = Clustering_moviemask.resultKMeans(iK).Vox_indCluster - matClusterID_vox;
        tempMat(tempMat~=0) = NaN;
        tempMat(~isnan(tempMat)) = 1;
        tempMat(isnan(tempMat))=0;
        matProb(:,iVox) = sum(tempMat, 2)./100;
    end
%     toc;

    matProb = single(matProb);
    save(fullfile(dirDataBOLD, sprintf('ClusteringProbability_TorRhoSigSpiArtMovie123_moviemask_%dMeans.mat',targetK)),...
        'Clustering_moviemask', 'matProb', '-v7.3')

    fprintf(1, 'K = %d: movie mask: Results saved \n', targetK);
% end
% 
% for iK=1:length(paramClustering_global.setK)
% 
%     targetK = paramClustering_global.setK(iK);
        fprintf(1, 'K = %d: brain mask \n', targetK);

%     tic;
    matProbBrain = NaN(size(Clustering_brainmask.resultKMeans(1).Vox_indCluster, 1));
    parfor iVox = 1:size(Clustering_brainmask.resultKMeans(1).Vox_indCluster, 1)
        matClusterID_vox = repmat(Clustering_brainmask.resultKMeans(iK).Vox_indCluster(iVox,:),...
            size(Clustering_brainmask.resultKMeans(iK).Vox_indCluster, 1), 1);
        tempMat = [];
        tempMat = Clustering_brainmask.resultKMeans(iK).Vox_indCluster - matClusterID_vox;
        tempMat(tempMat~=0) = NaN;
        tempMat(~isnan(tempMat)) = 1;
        tempMat(isnan(tempMat))=0;
        matProbBrain(:,iVox) = sum(tempMat, 2)./100;
    end
%     toc;
    
    matProbBrain = single(matProbBrain);
    save(fullfile(dirDataBOLD, sprintf('ClusteringProbability_TorRhoSigSpiArtMovie123_brainmask_%dMeans.mat',targetK)),...
        'Clustering_brainmask', 'matProbBrain', '-v7.3')
    fprintf(1, 'K = %d: brain mask: Results saved \n', targetK);

end

delete(pool);
