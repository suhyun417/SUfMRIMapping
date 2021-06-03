% doClusteringCorrMap_ROI_FaceCells.m
%
% 2021/05/17 SHP 
% Clustering of only face-selective cells based on ROI-correlation map

clear all;


dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

% Load correlation matrix of all cortical cells
load('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat', ...
    'Clustering_*')

% get only face selective cells
load('/procdata/parksh/_macaque/multipleFP_fsi.mat')
locFaceCell =  find(fsi.matFSI(:,1)>0.33); % find(abs(fsi.matFSI(:,1))>0.33);

meanCorrROI = Clustering_meanROI.matR(locFaceCell, :); % only for face selective neurons
absmaxCorrROI = Clustering_maxabsROI.matR(locFaceCell, :); 
catSubjID = Clustering_meanROI.catSubjID(locFaceCell, :); 
catAreaID = Clustering_meanROI.catAreaID(locFaceCell, :); 
catChanID = Clustering_meanROI.catChanID(locFaceCell, :);

nameROI = Clustering_meanROI.nameROI; 
setArea = Clustering_meanROI.setArea;

clear Clustering_*

%% Quick clustering on ROI map 
numRepeat = 100; % number of repetition for entire clustering

setK = 2:30; %15;
opts = statset('Display','final');
numReplicates = 5; %

paramClustering_faceCells.methods = 'KMeans';
paramClustering_faceCells.setK = setK;
paramClustering_faceCells.numReplicates = numReplicates;
paramClustering_faceCells.descriptions = 'Clustering of neurons based on the masked whole-brain correlation maps: multiple repetitions (numRepeat) of an execution of kmeans function, which had its own multiple "replicates" (numReplicates)';
paramClustering_faceCells.numRepeat = numRepeat;


flagSave = 1;

% % totalSS
[a, c, totalSS] = kmeans(meanCorrROI, 1); % single unit
Clustering_meanROI_faceCells.totalSS_SU = totalSS;
[a, c, totalSS] = kmeans(meanCorrROI', 1); % roi
Clustering_meanROI_faceCells.totalSS_roi = totalSS;

for iK = 1:length(setK)
    
    K = setK(iK);                
    
    SU_indCluster = NaN(size(meanCorrROI, 1), numRepeat);
    SU_sumD = NaN(K, numRepeat);
    
    roi_indCluster = NaN(size(meanCorrROI, 2), numRepeat);
    roi_sumD = NaN(K, numRepeat);

    for iRep = 1:numRepeat
        
        fprintf(1, ':: K = %d; mean corr for each ROI ::\n', K);
        
        % Cluster single units based on whole brain correlation
        [IDX_SU, C_SU, SUMD_SU] = kmeans(meanCorrROI, K,...
            'Replicates', numReplicates, 'Options', opts);
        
        SU_indCluster(:, iRep) = IDX_SU;
        SU_sumD(:, iRep) = SUMD_SU;
        
        % Cluster ROIs based on single unit correlation
        [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(meanCorrROI', K,...
            'Replicates', numReplicates, 'Options', opts);
        
        roi_indCluster(:, iRep) = IDX_Vox;
        roi_sumD(:, iRep) = SUMD_Vox;
        
    end

    Clustering_meanROI_faceCells.resultKMeans(iK).SU_indCluster = SU_indCluster;
    Clustering_meanROI_faceCells.resultKMeans(iK).SU_sumD = SU_sumD;
    Clustering_meanROI_faceCells.resultKMeans(iK).roi_indCluster = roi_indCluster;
    Clustering_meanROI_faceCells.resultKMeans(iK).roi_sumD = roi_sumD;   
    
%     
%     if flagSave
%         save(fullfile(dirDataBOLD, 'Clustering_CorrMap_Movie123_ArtRHROI_probability.mat'),...
%             'Clustering*', 'paramClustering*');
%         fprintf(1, ':: K = %d; Movie-driven mask :: Results saved \n', K);
%     end
end
Clustering_meanROI_faceCells.matR = meanCorrROI;
Clustering_meanROI_faceCells.catSubjID = catSubjID;
Clustering_meanROI_faceCells.catAreaID = catAreaID;
Clustering_meanROI_faceCells.catChanID = catChanID;
Clustering_meanROI_faceCells.setArea = setArea; 
Clustering_meanROI_faceCells.nameROI = nameROI;



% maximum (in absolute) corr for each ROI
[a, c, totalSS] = kmeans(absmaxCorrROI, 1); % single unit
Clustering_maxabsROI_faceCells.totalSS_SU = totalSS;
[a, c, totalSS] = kmeans(absmaxCorrROI', 1); % roi
Clustering_maxabsROI_faceCells.totalSS_roi = totalSS;
for iK = 1:length(setK)
    
    K = setK(iK);                
    
    % max of abs
    SU_indCluster = NaN(size(absmaxCorrROI, 1), numRepeat);
    SU_sumD = NaN(K, numRepeat);
    
    roi_indCluster = NaN(size(absmaxCorrROI, 2), numRepeat);
    roi_sumD = NaN(K, numRepeat);

    for iRep = 1:numRepeat
        
        fprintf(1, ':: K = %d; maxabs ::\n', K);
        
        % Cluster single units based on whole brain correlation
        [IDX_SU, C_SU, SUMD_SU] = kmeans(absmaxCorrROI, K,...
            'Replicates', numReplicates, 'Options', opts);
        
        SU_indCluster(:, iRep) = IDX_SU;
        SU_sumD(:, iRep) = SUMD_SU;
        
        % Cluster ROIs based on single unit correlation
        [IDX_Vox, C_Vox, SUMD_Vox] = kmeans(absmaxCorrROI', K,...
            'Replicates', numReplicates, 'Options', opts);
        
        roi_indCluster(:, iRep) = IDX_Vox;
        roi_sumD(:, iRep) = SUMD_Vox;
        
    end

    Clustering_maxabsROI_faceCells.resultKMeans(iK).SU_indCluster = SU_indCluster;
    Clustering_maxabsROI_faceCells.resultKMeans(iK).SU_sumD = SU_sumD;
    Clustering_maxabsROI_faceCells.resultKMeans(iK).roi_indCluster = roi_indCluster;
    Clustering_maxabsROI_faceCells.resultKMeans(iK).roi_sumD = roi_sumD;

end
Clustering_maxabsROI_faceCells.matR = absmaxCorrROI;
Clustering_maxabsROI_faceCells.catSubjID = catSubjID;
Clustering_maxabsROI_faceCells.catAreaID = catAreaID;
Clustering_maxabsROI_faceCells.catChanID = catChanID; 
Clustering_maxabsROI_faceCells.setArea = setArea; %{corrMap_Area(setArea).nameArea};
Clustering_maxabsROI_faceCells.nameROI = nameROI;
    
if flagSave
    save('/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat', ...
    '*_faceCells', '-append')
end

%% Explained variance elbow plot
setK = paramClustering_faceCells.setK; %Clustering.setK;

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
    
    matWSS(:,iK) = sum(Clustering_meanROI_faceCells.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
%     matExpVar(iK,1) = sum(tExpVar);

end

totalSS = Clustering_meanROI_faceCells.totalSS_SU;
% [a, c, totalSS] = kmeans(matR_SU_moviemask', 1); %kmeans(matR_SU', 1);
betweenSS = totalSS-matWSS;
% totalVar = totalSS/(2*size(matR_SU,2)); %totalD/(2*size(matR_SU,2));

propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

%% plot Clustering results
% explained variance
figure;
plot(setK, propExplained'.*100, 'ko-'); hold on
xlabel('Number of cluster (K)')
ylabel('Explained variance (%)')
set(gca, 'XTick', setK)

figure;
plot(setK(1:end-1), diff(propExplained').*100)
hold on
plot(setK(1:end-1), mean(diff(propExplained').*100, 2), 'ko-', 'LineWidth', 2)
title('difference of explained variance for each K: using mean ROI')

% distance within each cluster
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iK=1:length(setK)
    plot(iK+1, mean(Clustering_meanROI_faceCells.resultKMeans(iK).SU_sumD, 2), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
end
xlim([1 setK(end)])
set(gca, 'LineWidth', 2, 'FontSize', 12)
box off
title('Clustering of single units')
xlabel('Number of cluster')
ylabel('Within-cluster distance')