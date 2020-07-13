% doClusteringCorrMap_voxel_multipleFP_probability.m
%
% 2020/07/09 SHP
% Perform K-means clustering on correlation maps of 389 neurons from 4 face
% patches
% Repeat clustering procedure for each K 100 times so that we can have
% probability of clustering

clear all;

%% Settings
% Load correlation matrix of all cortical cells
nameSubjBOLD = 'Art';
load(sprintf('/procdata/parksh/_macaque/CorrMap_SU_AllCells%s_corticalFPMerged.mat', nameSubjBOLD), 'info*', 'corrMap_merged_FP'); 

% Load the masks (brain-only mask)
dirDataBOLD = sprintf('/procdata/parksh/_macaque/%s/', nameSubjBOLD);
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');

%
dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';
flagSave = 1;

% Only use voxels within brain
[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;
brainmask_vec = reshape(movieDrivenAmp.map_sm_brain>0, nVox, 1); % change the 3D mask to 1D

matR_SU_all_brainmask = corrMap_merged_FP.matR(brainmask_vec,:); % 27113 voxels; 
catSubjID = corrMap_merged_FP.catSubjID;
catAreaID = corrMap_merged_FP.catAreaID;


%% Clustering of voxles based on unthresholded map 
numRepeat = 100; % number of repetition for entire clustering

setK = 2:40; %15;
opts = statset('Display','final');
numReplicates = 5; %

paramClustering_global.methods = 'KMeans';
paramClustering_global.setK = setK;
paramClustering_global.numReplicates = numReplicates;
paramClustering_global.descriptions = 'Clustering of neurons based on the masked whole-brain correlation maps: multiple repetitions (numRepeat) of an execution of kmeans function, which had its own multiple "replicates" (numReplicates)';
paramClustering_global.numRepeat = numRepeat;


% % totalSS
[a, c, totalSS] = kmeans(matR_SU_all_brainmask, 1);
Clustering_brainmask.totalSS = totalSS;


for iK = 1:length(setK)
    
    K = setK(iK);                
    
    Vox_indCluster_brainmask = NaN(sum(brainmask_vec), numRepeat);
    Vox_sumD_brainmask = NaN(K, numRepeat);

    for iRep = 1:numRepeat
        
        fprintf(1, ':: K = %d; voxel clustering ::\n', K);
        
        % Cluster voxels based on 389 singel unit correlation
        [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(matR_SU_all_brainmask, K, 'Replicates', numReplicates, 'Options', opts); %, numReplicates, 'Options', opts);
        
        Vox_indCluster_brainmask(:, iRep) = IDX_Vox;
        Vox_sumD_brainmask(:, iRep) = SUMD_Vox;
        
    end

    Clustering_brainmask.resultKMeans(iK).Vox_indCluster = Vox_indCluster_brainmask;
    Clustering_brainmask.resultKMeans(iK).Vox_sumD = Vox_sumD_brainmask;
    
    if flagSave
        save(fullfile(dirDataBOLD, sprintf('Clustering_4FPCells%sMovie123_Voxel_probability.mat', nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
            'Clustering*', 'paramClustering*');
         fprintf(1, ':: K = %d; voxel clustering: brainmask :: Results saved \n', K);
    end

end
paramClustering_brainmask.matR = matR_SU_all_brainmask;
paramClustering_brainmask.catSubjID = catSubjID;
paramClustering_brainmask.catAreaID = catAreaID;
paramClustering_brainmask.setArea = corrMap_merged_FP.setArea;

if flagSave
    save(fullfile(dirDataBOLD, sprintf('Clustering_4FPCells%sMovie123_Voxel_probability.mat', nameSubjBOLD)), ... 
       'paramClustering_brainmask', '-append');
    fprintf(1, ':: Voxel clustering: brainmask is DONE! :: Results saved in %s \n', ...
        fullfile(dirDataBOLD, sprintf('Clustering_4FPCells%sMovie123_Voxel_probability.mat', nameSubjBOLD)));
end

% 
% 
% %% Explained variance elbow plot
% setK = paramClustering_global.setK; %Clustering.setK;
% 
% matWSS=[];
% matExpVar=[];
% for iK = 1:length(setK)
%     curK = setK(iK);
% %     indClust = Clustering_moviemask_valid.resultKMeans(iK).SU_indCluster; %Clustering.resultKMeans(iK).SU_indCluster;
% %     [sortedClust, indSortedChan]=sort(indClust);
%     
% %     tExpVar=[];
% %     for ii = 1:curK
% %         tExpVar(ii,1) = Clustering.resultKMeans(iK).SU_sumD(ii)/(2*sum(sortedClust==ii));
% %     end
%     
%     matWSS(:,iK) = sum(Clustering_brainmask.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
% %     matExpVar(iK,1) = sum(tExpVar);
% 
% end
% 
% totalSS = Clustering_brainmask.totalSS_SU;
% % [a, c, totalSS] = kmeans(matR_SU_moviemask', 1); %kmeans(matR_SU', 1);
% betweenSS = totalSS-matWSS;
% % totalVar = totalSS/(2*size(matR_SU,2)); %totalD/(2*size(matR_SU,2));
% 
% propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;
% 
% %% plot Clustering results
% % explained variance
% figure;
% plot(setK, propExplained'.*100, 'ko-'); hold on
% xlabel('Number of cluster (K)')
% ylabel('Explained variance (%)')
% set(gca, 'XTick', setK)
% 
% figure;
% plot(setK(1:end-1), diff(propExplained').*100)
% hold on
% plot(setK(1:end-1), mean(diff(propExplained').*100, 2), 'ko-', 'LineWidth', 2)
% title('difference of explained variance for each K: using mean ROI')
% 
% % distance within each cluster
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% for iK=1:length(setK)
%     plot(iK+1, mean(Clustering_brainmask.resultKMeans(iK).SU_sumD, 2), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
%     hold on
% end
% xlim([1 setK(end)])
% set(gca, 'LineWidth', 2, 'FontSize', 12)
% box off
% title('Clustering of single units')
% xlabel('Number of cluster')
% ylabel('Within-cluster distance')
% 

%%





% %% Perform gap statistics to get the optimal number of clusters
% % Gap statistics to get the optimal number
% fprintf(1, 'Performing gap statistics on SU clustering... \n');
% eva_clusteringSU = evalclusters(matR_roi, 'kmeans', 'gap', 'KList', paramClustering_global.setK, 'Distance', 'sqEuclidean');
% save(fullfile(dirDataBOLD, 'Clustering_CorrMap_Movie123_ArtRHROI_probability_eval.mat'), 'eva_clusteringSU');
% fprintf(1, '...done! Evaluation results saved. \n')
% 
% % fprintf(1, 'Performing gap statistics on Voxel clustering... \n');
% % eva_clusteringVox = evalclusters(matR_roi', 'kmeans', 'gap', 'KList', paramClustering_global.setK, 'Distance', 'sqEuclidean'); 
% % fprintf(1, '...done! Evaluation results saved. \n')
% % save(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1_eval.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)), ... %nameSubjNeural, nameSubjBOLD)),...
% %         'eva_clusteringSU', 'eva_clusteringVox');
% 
% if flagParallel
%     delete(pool);
% end







