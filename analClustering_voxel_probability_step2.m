% analClustering_voxel_probability_step2.m

% 2017/01/29 SHP
% 1) For each K, load a matrix of likelihood that is computed by
% "analClustering_voxel_probability.m"
% 2) Visualize the results

clear all;

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjBOLD = 'Art';
dirDataHome = '/procdata/parksh/';
dirDataBOLD= fullfile(dirDataHome, nameSubjBOLD); %fullfile('/procdata/parksh/', nameSubjBOLD); %

% load voxel clustering results
load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  
% load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_voxel_probability_critCorr2.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  
 
%% Compute the proportion of total variance explained
setK = paramClustering_global.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
%     indClust = Clustering_moviemask_valid.resultKMeans(iK).SU_indCluster; %Clustering.resultKMeans(iK).SU_indCluster;
%     [sortedClust, indSortedChan]=sort(indClust);
    
%     tExpVar=[];
%     for ii = 1:curK
%         tExpVar(ii,1) = Clustering.resultKMeans(iK).SU_sumD(ii)/(2*sum(sortedClust==ii));
%     end
    
    matWSS(:,iK) = sum(Clustering_moviemask_valid.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
%     matExpVar(iK,1) = sum(tExpVar);

end

totalSS = Clustering_moviemask_valid.totalSS_SU;
% [a, c, totalSS] = kmeans(matR_SU_moviemask', 1); %kmeans(matR_SU', 1);
betweenSS = totalSS-matWSS;
% totalVar = totalSS/(2*size(matR_SU,2)); %totalD/(2*size(matR_SU,2));
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

%% 
K=7;
locMode = find(propExplained(:,K-1)==mode(propExplained(:,K-1)));
locMin = find(propExplained(:,K-1)==min(propExplained(:,K-1)));
[sortClust_validVoxel sortCell_validVoxel] = sort(Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode(1)));
figure;
imagesc(Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(sortCell_validVoxel,locMode)) % check whether clustering is identical in all those cases

figure;
imagesc(Clustering_moviemask_valid.resultKMeans(K-1).matProb(sortCell_validVoxel, sortCell_validVoxel)); % how probable


% clear all;
% 
% 
% % sort voxels based on 2-means clustering
% K=2;
% load(sprintf('/procdata/parksh/Art/ClusteringProbability_TorRhoSigSpiArtMovie123_moviemask_%dMeans.mat', K))
% 
% [a, sortedIndVox] = sort(matProb(:,1), 'descend');
% 
% K=4;
% load(sprintf('/procdata/parksh/Art/ClusteringProbability_TorRhoSigSpiArtMovie123_moviemask_%dMeans.mat', K))
% figure;
% imagesc(matProb(sortedIndVox, sortedIndVox))
% 
% 
% % Check each K for how stable  the "best" clustering that computed by K-means was
% figCheckK = figure;
% set(figCheckK, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 880 765])
% 
% % load the Clustering structure for just once
% load('/procdata/parksh/Art/ClusteringProbability_TorRhoSigSpiArtMovie123_moviemask_2Means.mat', 'Clustering_moviemask')
% 
% for K=2:20
% 
%     load(sprintf('/procdata/parksh/Art/ClusteringProbability_TorRhoSigSpiArtMovie123_moviemask_%dMeans.mat', K), 'matProb')
% %     imagesc(matProb(sortedIndVox, sortedIndVox))
% %     [aa, bb] = min(sum(Clustering_moviemask.resultKMeans(K-1).Vox_sumD)); % the best answer by K-means
% %     [sortedClust, sortedIndVox_cluster] = sort(Clustering_moviemask.resultKMeans(K-1).Vox_indCluster(:,bb), 'descend');
%     
%     figure(figCheckK); clf;
%     imagesc(matProb) %imagesc(matProb(sortedIndVox_cluster, sortedIndVox_cluster));
%     set(gca, 'CLim', [0 1])
%     colorbar;
%     title(sprintf('Likelihood of voxels being clustered together: K = %d', K))
%     
%     input('')
%     
% end
% 
% % % Let's use K=10 case as a reference
% % K=10;
% % load(sprintf('/procdata/parksh/Art/ClusteringProbability_TorRhoSigSpiArtMovie123_moviemask_%dMeans.mat', K), 'matProb')
% % %     imagesc(matProb(sortedIndVox, sortedIndVox))
% %     [aa, bb] = min(sum(Clustering_moviemask.resultKMeans(K-1).Vox_sumD)); % the best answer by K-means
% %     [sortedClust, sortedIndVox_cluster] = sort(Clustering_moviemask.resultKMeans(K-1).Vox_indCluster(:,bb), 'descend');
% %     
% % for K=2:20
% % 
% %     load(sprintf('/procdata/parksh/Art/ClusteringProbability_TorRhoSigSpiArtMovie123_moviemask_%dMeans.mat', K), 'matProb')
% % %     imagesc(matProb(sortedIndVox, sortedIndVox))
% % %     [aa, bb] = min(sum(Clustering_moviemask.resultKMeans(K-1).Vox_sumD)); % the best answer by K-means
% % %     [sortedClust, sortedIndVox_cluster] = sort(Clustering_moviemask.resultKMeans(K-1).Vox_indCluster(:,bb), 'descend');
% %     
% %     figure(figCheckK); clf;
% %     imagesc(matProb(sortedIndVox_cluster, sortedIndVox_cluster));
% %     set(gca, 'CLim', [0 1])
% %     colorbar;
% %     title(sprintf('Likelihood of voxels being clustered together: K = %d', K))
% %     
% %     input('')
% %     
% % end
    


