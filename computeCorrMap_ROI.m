% computeCorrMap_ROI.m
%
% 2020/06/29

clear all;

% Load the ROI indices
nameSubjBOLD = 'Art';
load(sprintf('/procdata/parksh/_macaque/%s/ROIs/%s_ROIs_set00_RH.mat', nameSubjBOLD));

% Load correlation matrix of all cortical cells
load(sprintf('/procdata/parksh/_macaque/CorrMap_SU_AllCells%s_corticalFPMerged.mat', nameSubjBOLD), 'info*', 'corrMap_Area'); %, 'info*', 'corrMap_Area', 'corrMap_merged');

dirFig = '/projects/parksh/NeuroMRI/_labNote/_figs';

numROI = size(matROIIndices, 2);
nVoxPerROI = sum(matROIIndices);

setArea = 1:4 ;%1:5;
matR_SU_all = cat(2, corrMap_Area(setArea).matR); % except for no face patch neurons
catSubjID = cat(1, corrMap_Area(setArea).catSubjID);
catAreaID = floor(cat(1, corrMap_Area(setArea).catSubjID)./10);

sumCorrROI = matR_SU_all'*matROIIndices;
meanCorrROI = sumCorrROI./repmat(nVoxPerROI, size(sumCorrROI, 1), 1);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1700 390])
imagesc(meanCorrROI')
set(gca, 'XTick', cat(1, find(diff(catAreaID)>0), size(meanCorrROI, 1)), 'XTickLabel', {corrMap_Area(setArea).nameArea})
set(gca, 'YTick', 1:numROI, 'YTickLabel', paramROI.nameROI)
set(gca, 'CLim', [-1 1].*0.5)
set(gca, 'TickDir', 'out', 'Box', 'off')
colorbar;
title('mean correlation value for each ROI')

for iRoi = 1:numROI
    setCorr = [];
    setCorr = matR_SU_all(matROIIndices(:, iRoi)>0, :);
    [y, i] = max(abs(setCorr));
    ind = sub2ind(size(setCorr), i, 1:size(setCorr,2));
    
    corrROI(iRoi).setCorr = setCorr;
    corrROI(iRoi).mean = mean(setCorr);
    corrROI(iRoi).max = max(setCorr);
    corrROI(iRoi).absmax = setCorr(ind);
    
end

maxCorrROI = cat(1, corrROI.max)';
absmaxCorrROI = cat(1, corrROI.absmax)';

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1700 390])
imagesc(absmaxCorrROI')
set(gca, 'XTick', cat(1, find(diff(catAreaID)>0), size(meanCorrROI, 1)), 'XTickLabel', {corrMap_Area(setArea).nameArea})
set(gca, 'YTick', 1:numROI, 'YTickLabel', paramROI.nameROI)
set(gca, 'CLim', [-1 1].*0.8)
set(gca, 'TickDir', 'out', 'Box', 'off')
colorbar;
title('max correlation value for each ROI')


figROI = figure;
set(figROI, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position',  [100 500 785 408])
for iR = 1:numROI
    figure(figROI);
    clf;
    hist(corrROI(iR).setCorr(:, randperm(389, 10)), 20);
    xlim([-1 1].*0.6)
    title(sprintf('ROI %d: %s', iR, paramROI.nameROI{iR}));
    input('')
end

%% Quick clustering on masked map 
dirDataBOLD = '/procdata/parksh/_macaque/Art/';

numRepeat = 100; % number of repetition for entire clustering

setK = 2:20; %15;
opts = statset('Display','final');
numReplicates = 5; %

paramClustering_global.methods = 'KMeans';
paramClustering_global.setK = setK;
paramClustering_global.numReplicates = numReplicates;
paramClustering_global.descriptions = 'Clustering of neurons based on the masked whole-brain correlation maps: multiple repetitions (numRepeat) of an execution of kmeans function, which had its own multiple "replicates" (numReplicates)';
paramClustering_global.numRepeat = numRepeat;

% if flagParallel
%     pool = parpool;                      % Invokes workers
%     stream = RandStream('mlfg6331_64');  % Random number stream
%     opts = statset('UseParallel',1,'UseSubstreams',1, 'Streams',stream,...
%         'MaxIter', 1000, 'Display','final');
%     paramClustering_global.parallel = opts;
% end

% [nx ny nz] = size(movieDrivenAmp.mask_amp1);
% nVox = nx*ny*nz;


flagSave = 1;
% % totalSS
[a, c, totalSS] = kmeans(meanCorrROI, 1); % single unit
Clustering_meanROI.totalSS_SU = totalSS;
[a, c, totalSS] = kmeans(meanCorrROI', 1); % roi
Clustering_meanROI.totalSS_roi = totalSS;


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

    Clustering_meanROI.resultKMeans(iK).SU_indCluster = SU_indCluster;
    Clustering_meanROI.resultKMeans(iK).SU_sumD = SU_sumD;
    Clustering_meanROI.resultKMeans(iK).roi_indCluster = roi_indCluster;
    Clustering_meanROI.resultKMeans(iK).roi_sumD = roi_sumD;   
    
%     
%     if flagSave
%         save(fullfile(dirDataBOLD, 'Clustering_CorrMap_Movie123_ArtRHROI_probability.mat'),...
%             'Clustering*', 'paramClustering*');
%         fprintf(1, ':: K = %d; Movie-driven mask :: Results saved \n', K);
%     end
end
Clustering_meanROI.matR = meanCorrROI;
Clustering_meanROI.catSubjID = catSubjID;
Clustering_meanROI.catAreaID = catAreaID;
Clustering_meanROI.setArea = {corrMap_Area(setArea).nameArea};
Clustering_meanROI.nameROI = paramROI.nameROI;



% maximum (in absolute) corr for each ROI
[a, c, totalSS] = kmeans(absmaxCorrROI, 1); % single unit
Clustering_maxabsROI.totalSS_SU = totalSS;
[a, c, totalSS] = kmeans(absmaxCorrROI', 1); % roi
Clustering_maxabsROI.totalSS_roi = totalSS;
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

    Clustering_maxabsROI.resultKMeans(iK).SU_indCluster = SU_indCluster;
    Clustering_maxabsROI.resultKMeans(iK).SU_sumD = SU_sumD;
    Clustering_maxabsROI.resultKMeans(iK).roi_indCluster = roi_indCluster;
    Clustering_maxabsROI.resultKMeans(iK).roi_sumD = roi_sumD;

end
Clustering_maxabsROI.matR = absmaxCorrROI;
Clustering_maxabsROI.catSubjID = catSubjID;
Clustering_maxabsROI.catAreaID = catAreaID;
Clustering_maxabsROI.setArea = {corrMap_Area(setArea).nameArea};
Clustering_maxabsROI.nameROI = paramROI.nameROI;
    
if flagSave
    save(fullfile(dirDataBOLD, 'Clustering_CorrMap_4FPs_Movie123_ArtRHROI_probability.mat'),...
        'Clustering*', 'paramClustering*');
%     fprintf(1, ':: K = %d; Movie-driven mask :: Results saved \n', K);
end

%% Explained variance elbow plot
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
    
    matWSS(:,iK) = sum(Clustering_meanROI.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
%     matExpVar(iK,1) = sum(tExpVar);

end

totalSS = Clustering_meanROI.totalSS_SU;
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
    plot(iK+1, mean(Clustering_meanROI.resultKMeans(iK).SU_sumD, 2), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
end
xlim([1 setK(end)])
set(gca, 'LineWidth', 2, 'FontSize', 12)
box off
title('Clustering of single units')
xlabel('Number of cluster')
ylabel('Within-cluster distance')


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







