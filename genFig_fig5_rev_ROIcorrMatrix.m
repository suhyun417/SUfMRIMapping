% genFig_fig5_rev.m
%
% 2017/02/21 SHP
% Revised figure 5 to show ROI-based analysis of correlation maps of each
% cell group

clear all;

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
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');
    
% % Add necessary toolbox 
% addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution

% Set directories 
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjNeural = 'Spi'; % 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Directory for saving figures as graphic files
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');

%% Load data
% 1) Clustering results
load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  
% load(fullfile(dirDataNeural, sprintf('Clustering_%s%sMovie123_new_masked.mat', nameSubjNeural, nameSubjBOLD)));
% 2) Movie-driven mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp');
% 3) Valid cells considering likelihood being clustered together
load(fullfile(dirDataNeural, 'Clustering_SU_allCells_validVoxels_critCorr1_7Means_prob.mat'))


%% Load the corr map and select valid channel: subject by subject
numSubject = size(setNameSubjNeural, 2);
matR_SU_all = [];

for iSubj = 1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');
    
    % 2017/01/19: Excluded some units, in case of Spice's data which hasn't been went through this procedure yet
    switch lower(nameSubjNeural)
        case 'spi'
            excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
            validChanIndex_clustering = setdiff(paramCorr.validChanIndex, excChanIndex);
            validChanID_clustering = cat(2, paramCorr.validChanID(validChanIndex_clustering,:), num2str(ones(size(validChanIndex_clustering)).*iSubj)); %paramCorr.validChanID(validChanIndex_clustering,:);
        otherwise
            validChanIndex_clustering = (1:length(paramCorr.validChanIndex))';
            validChanID_clustering = cat(2, paramCorr.validChanID, num2str(ones(size(paramCorr.validChanIndex)).*iSubj)); %paramCorr.validChanID;
    end
    
    matR_SU_valid = matR_SU(:, validChanIndex_clustering);
    matR_SU_all = cat(2, matR_SU_all, matR_SU_valid);
    
     clear matR_SU matR_SU_valid
end

[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;

% Apply movie-driven mask to correlation matrix 
% 2017/01/19 & valid channels 
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D
matR_SU_all_moviemask = matR_SU_all(moviemask_vec,:); %matR_SU(moviemask_vec,:); % 15495 voxels
matR_SU_all_moviemask_valid  = matR_SU_all_moviemask(paramClustering_global.locValidVox, :);

%% Get the ROIs and compute the average correlation for each cell
% From a manually drawn ROI + face patches
typeMotionROI = 2; %1; % type 2 contains more restricted motion-area-ROI (within STS) than type 1
load(sprintf('/procdata/parksh/Art/ROIs/Art_ROIs_RetinoMotion%dFace_bothH.mat', typeMotionROI), 'matROIIndices', 'paramROI')

matROIIndices_moviemask = matROIIndices(moviemask_vec, :);
matROIIndices_moviemask_valid  = matROIIndices_moviemask(paramClustering_global.locValidVox, :);

nVoxPerROI = sum(matROIIndices_moviemask_valid);

sumCorrROI = matR_SU_all_moviemask_valid'*matROIIndices_moviemask_valid;
meanCorrROI = sumCorrROI./repmat(nVoxPerROI, size(sumCorrROI, 1), 1);


%% Compute the proportion of total variance explained to get the clustering index
setK = paramClustering_global.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)    
    matWSS(:,iK) = sum(Clustering_moviemask_valid.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
end

totalSS = Clustering_moviemask_valid.totalSS_SU;
betweenSS = totalSS-matWSS;

propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

%% FIg 5A: Correlation matrix between cells and ROIs
K=7;
locMode = find(propExplained(:,K-1)==mode(propExplained(:,K-1)));
indClust = Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode(1));
% locMin = find(propExplained(:,K-1)==min(propExplained(:,K-1)));
% [sortClust_validVoxel, sortCell_validVoxel] = sort(Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode(1)));

orderROI = [1 3 2]; % paramROI.nameROI 1. V1+ / 2. Motion areas / 3. Face patches  %[4 2 5]; %[4 1 3 2 5]; %% 1: parafoveal:ventral stream / 2: face patches / 3: Eccentric visual cortex / 4: Central V1 / 5: MT
orderClust_Cell = [1 2 5 6 3 7 4]; %[4 2 5 1 7 6 3]; %1:numClust_Cell; %[3 2 7 6 1 5 4]; %[7 6 2 3 4 8 1 5]; %[4 1 2 5 6 7 3]; %  [6 5 4 3 2 1 7 ]; %[6 5 3 2 4 1 7 ]; % from pos to neg
orderClustValid = [4 2 3 6 7 5 1]; 

catIndValidCells = cat(1, resultProbClustering(orderClustValid).validIndCells);
% reOrdSortCell=[];
% for iCC=1:K
%     reOrdSortCell = cat(1, reOrdSortCell, find(indClust==orderClust_Cell(iCC)));
% end

 figCorrMat_voxCluster = figure;
 set(figCorrMat_voxCluster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1550 330])
 imagesc(meanCorrROI(catIndValidCells, orderROI)')
 set(gca, 'CLim', [-1 1].*0.5)
%  colorbar('EastOutside');
 set(gca, 'YTick', 1:length(orderROI), 'YTickLabel', paramROI.nameROI(orderROI))
 set(gca, 'XTick', 1:135, 'XTickLabel', num2str(indClust(catIndValidCells)))
 set(gca, 'TickDir', 'out')
 box off
 
 print(figCorrMat_voxCluster,...
     fullfile(dirFig, 'figRev5A_orderClusterSize'),'-dtiff', '-r200')


%% Fig 5B: average correlation between cell groups and ROIs
K=7;
locMode = find(propExplained(:,K-1)==mode(propExplained(:,K-1)));
indClust = Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode(1));

reOrdMeanCorrROI=[];
orderClustValid = [4 2 3 6 7 5 1];
for iK=1:K    
    iCC = orderClustValid(iK);
    reOrdMeanCorrROI = cat(1, reOrdMeanCorrROI, meanCorrROI(resultProbClustering(iCC).validIndCells, :));
    
    meanCorrROI_Cluster(iK,:) = mean(meanCorrROI(resultProbClustering(iCC).validIndCells, :));
    steCorrROI_Cluster(iK,:) = std(meanCorrROI(resultProbClustering(iCC).validIndCells, :))./sqrt(length(resultProbClustering(iCC).validIndCells));
end

colorOrder3 = [142 50 233; 228 153 203; 65 225 150]./255;
figBarCG = figure;
set(figBarCG, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [200 200 1175 290])
BH = bar(meanCorrROI_Cluster(:,orderROI)); hold on;
set(BH, {'EdgeColor'}, mat2cell(colorOrder3, [1 1 1], [3]))
set(BH, 'FaceColor', 'none')
set(BH, 'LineWidth', 2)

xCoords=[];
for iBH = 1:length(BH)
    tempX = get(get(BH(iBH), 'Children'), 'XData'); % get the x coordinates of each bar
    xCoords(:, iBH) = mean(tempX(2:3,:));
end

figure(figBarCG);
hold on;
for iK = 1:7
    iCC = orderClustValid(iK);
    PH = plot(repmat(xCoords(iK, :), length(resultProbClustering(iCC).validIndCells), 1), meanCorrROI(resultProbClustering(iCC).validIndCells, orderROI),...
        'o', 'LineWidth', 1.5);
    set(PH, {'MarkerEdgeColor'}, mat2cell(colorOrder3, [1 1 1], [3]))
%     set(PH, 'MarkerFaceColor', 'w') %
    set(PH, {'MarkerFaceColor'}, mat2cell(colorOrder3, [1 1 1], [3]))
    set(PH, 'MarkerSize', 4)
end

set(gca, 'TickDir', 'out')
box off
ylim([-1 1].*0.6)
set(gca, 'XTick', 1:7, 'XTickLabel', orderClustValid)

 print(figBarCG, fullfile(dirFig, 'figRev5B_avgCorrROI_indDataPoint_orderClusterSize'), '-depsc')



% 
% %% Probability of the clustering
% K=7;
% locMode = find(propExplained(:,K-1)==mode(propExplained(:,K-1)));
% locMin = find(propExplained(:,K-1)==min(propExplained(:,K-1)));
% indClust = Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode(1)); %
% [sortClust_validVoxel, sortCell_validVoxel] = sort(Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode(1)));
% 
% orderROI = [4 2 5]; %[4 1 3 2 5]; %% 1: parafoveal:ventral stream / 2: face patches / 3: Eccentric visual cortex / 4: Central V1 / 5: MT
% orderClust_Cell = [4 2 5 1 7 6 3]; %
% 
% reOrdSortCell=[];
% for iCC=1:K
%     reOrdSortCell = cat(1, reOrdSortCell, find(indClust==orderClust_Cell(iCC)));
% end
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto'); %, 'Position', [100 100 200 950])
% imagesc(Clustering_moviemask_valid.resultKMeans(K-1).matProb(reOrdSortCell, reOrdSortCell)); % how probable


% %% Select out cells that are not stable
% K=7;
% matProb = Clustering_moviemask_valid.resultKMeans(K-1).matProb;
% indClust = Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode(1)); %
% for iCC=1:K
%     indCellCurClust = find(indClust==orderClust_Cell(iCC));
%     freqLowProb = sum(matProb(indCellCurClust, indCellCurClust)<0.5); % 50% of the clustering
%     locCells = find(freqLowProb < (length(indCellCurClust)/2)); % doesn't cluster with half of the cells in the current cluster
%     resultProbClustering(iCC).clusterID = iCC;
%     resultProbClustering(iCC).clusterID_org = orderClust_Cell(iCC);
%     resultProbClustering(iCC).validIndCells = indCellCurClust(locCells);
% end



% for iCC=1:K
%     figure;
%     set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [500 500 250 150])
%     b(iCC) = bar(1:3, meanCorrROI_Cluster(iCC, orderROI));
%     set(b(iCC), 


% % k-means clustering on cells to order the matrix
% eva = evalclusters(meanCorrROI, 'kmeans', 'gap', 'KList', [2:15], 'Distance', 'sqEuclidean');
% numClust_Cell = eva.OptimalK; % 7, according to gap statistics

% [indClust, C] = kmeans(meanCorrROI, numClust_Cell,  'Replicates', 5, 'Display', 'final');
% [Y, sortCell] = sort(indClust);

% numRepeat = 100;
% numClust_Cell=8;
% 
% SU_indCluster = NaN(size(matR_SU_all_moviemask_valid,2), numRepeat);
% SU_sumD = NaN(numClust_Cell, numRepeat);
% for iRep = 1:numRepeat    
%     % Cluster single units based on whole brain correlation
%     [IDX_SU, C_SU, SUMD_SU] = kmeans(meanCorrROI, numClust_Cell,...
%         'Replicates', 5,  'Display', 'final');
%     
%     SU_indCluster(:, iRep) = IDX_SU;
%     SU_sumD(:, iRep) = SUMD_SU;    
% end

% [aa, locBest] = min(sum(SU_sumD)); % the best answer by K-means
% [sortedClust, sortedCells_ROI] = sort(SU_indCluster(:,locBest));
% figure;
% imagesc(matProb(sortedCells_ROI, sortedCells_ROI))
% title(sprintf('likelihood of clustering: K=%d, using %d ROIs', numClust_Cell, size(meanCorrROI,2)))
% figure;
% imagesc(meanCorrROI(sortedCells_ROI,:))
% title(sprintf('Clustering of single units: K=%d, using %d ROIs', numClust_Cell, size(meanCorrROI,2)))

% orderROI = [4 1 3 2 5]; %[3 1 4 2]; %[3 4 5 1 2]; %[1 5 4 3 2]; % 1: face patch / 2: screen boundary? / 3: V1F / 4: V1PF / 5: MT
% orderClust_Cell = [6 3 4 8 7 5 1 2]; %1:numClust_Cell; %[3 2 7 6 1 5 4]; %[7 6 2 3 4 8 1 5]; %[4 1 2 5 6 7 3]; %  [6 5 4 3 2 1 7 ]; %[6 5 3 2 4 1 7 ]; % from pos to neg
% 
% indClust = SU_indCluster(:,locBest);
% reOrdMeanCorrROI=[];
% for iCC=1:numClust_Cell
%     reOrdMeanCorrROI = cat(1, reOrdMeanCorrROI, meanCorrROI(indClust==orderClust_Cell(iCC), :));
% end
% reOrdSortCell=[];
% for iCC=1:numClust_Cell
%     reOrdSortCell = cat(1, reOrdSortCell, find(indClust==orderClust_Cell(iCC)));
% end
% 
% figCorrMat_voxCluster = figure;
% set(figCorrMat_voxCluster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1550 330])
% imagesc(reOrdMeanCorrROI(:, orderROI)')
% set(gca, 'CLim', [-1 1].*0.5)
% colorbar('EastOutside');
% set(gca, 'YTick', 1:5, 'YTickLabel', orderROI)
% set(gca, 'XTick', 1:135, 'XTickLabel', num2str(indClust(reOrdSortCell)))
% 
% print(figCorrMat_voxCluster,...
%     fullfile(dirFig, sprintf('corrMat_allCells_5ROIs_validVox_critCorr1_%dMeans_V1FSort_hor', max(orderClust_Cell))),...
%     '-dtiff', '-r200')
% 
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 200 950])
% imagesc(reOrdMeanCorrROI(:, orderROI))
% set(gca, 'CLim', [-1 1].*0.5)
% colorbar('SouthOutside');
% set(gca, 'XTick', 1:5, 'XTicklabel', orderROI)
% set(gca, 'YTick', 1:135, 'YTickLabel', num2str(indClust(reOrdSortCell)))
% print(gcf, fullfile(dirFig, sprintf('corrMat_allCells_5ROIs_validVox_critCorr1_%dMeans_V1FSort_ver', max(orderClust_Cell))),...
%     '-dtiff', '-r200')

% % make blue-white-red colorbar
% cval = 0.5;
% cmin = -cval; cmax = cval;
% colornum = 256;
% colorInput = [1 0 0; 1 1 1; 0 0 1];
% oldSteps = linspace(-1, 1, length(colorInput));
% newSteps = linspace(-1, 1, colornum);
% for j=1:3 % RGB
%     newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
% end
% endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
% newmap = squeeze(newmap_all(1:endPoint, :));
% figure(gcf)
% set(gca, 'CLim', [cmin cmax])
% colormap(flipud(newmap))


for targetK = [6 8 10 12 13];
    
    nVoxPerROI=[]; sumCorrROI=[]; meanCorrROI=[];
    
    % Prep matrix of clustering indices for computation
    [aa, locBest] = min(sum(Clustering_moviemask_valid.resultKMeans(targetK-1).Vox_sumD)); % the best answer by K-means    

    matClusterInd = NaN(size(Clustering_moviemask_valid.resultKMeans(1).Vox_indCluster, 1), targetK);
    for iROI = 1:targetK
        matClusterInd(:,iROI) = (Clustering_moviemask_valid.resultKMeans(targetK-1).Vox_indCluster(:,locBest) == iROI);
    end
    
    nVoxPerROI = sum(matClusterInd);
    
    sumCorrROI = matR_SU_all_moviemask_valid'*matClusterInd;
    meanCorrROI = sumCorrROI./repmat(nVoxPerROI, size(sumCorrROI, 1), 1);
    
    % interim k-means to order cells
    [indClust, C] = kmeans(meanCorrROI, 5,  'Replicates', 5, 'Display', 'final');
    [Y, sortCell] = sort(indClust);
    
    figCorrMat_voxCluster = figure;
    set(figCorrMat_voxCluster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1165 175])
    imagesc(meanCorrROI')
    set(gca, 'CLim', [-1 1].*0.6)
    colorbar;
    title(sprintf('Correlation between neurons and voxels: %d ROIs', targetK))
    ylabel('Voxel clusters')
    xlabel('Cells')
    print(figCorrMat_voxCluster, fullfile(dirFig, sprintf('CorrMatrix_%s%s_%dROI', cell2mat(setNameSubjNeural), nameSubjBOLD, targetK)),...
        '-dtiff', '-r200')
    
end

setTargetK = [6 8 10 12 13];
clusterOrder_arbitrary = cell(size(setTargetK));
clusterOrder_arbitrary{1} = [1 2 4 6 3 5]; % K=6: 1: face patches / 2: V1+alpha / 4: caudal STS (MT, FST,..) + LIP + arcuate / 6: boundary of screen? / 3: medial frontal / 5: medial
clusterOrder_arbitrary{2} = [3 8 4 1 6 2 5 7]; % K=8: same first six as above + cluster 5 & 7 are junk (mainly divided from 3 and 5 from K=6 case)
clusterOrder_arbitrary{3} = [1 6 5 4 7 3 8 9 10 2]; % K=10: 1: face patches / 6: V1 / 5: V1PF? / 4: caudal STS (MT, FST...) + LIP + arcuate / 7: boundary of screen? / [3 8 9 10 2] junk
clusterOrder_arbitrary{4} = [8 4 3 9 7 10 12 11 6 5 1 2]; % K=12; 8: face patches / 4: V1 / 3: V1PF? / 9: MT, FST / 7: ? divided from caudal STS chunk / 10: boundary of screen? / [12 11 6 5 1 2] junk
clusterOrder_arbitrary{5} = [1 11 4 12 6 8 10 7 9 2 5 13 3]; % K=13; 1: face patches / 11: V1 / 4: V1PF? / 12: MT, FST / 6: ? divided from caudal STS chunk / 8: boundary of screen? / [10 7 9 2 5 13 3] junk


for iK = 1:length(setTargetK);
    targetK = setTargetK(iK);
    nVoxPerROI=[]; sumCorrROI=[]; meanCorrROI=[];
    
    % Prep matrix of clustering indices for computation
    matClusterInd = NaN(length(Clustering_moviemask.resultKMeans(1).Vox_indCluster), targetK);
    for iROI = 1:targetK
        matClusterInd(:,iROI) = (Clustering_moviemask.resultKMeans(targetK-1).Vox_indCluster == iROI);
    end
    
    nVoxPerROI = sum(matClusterInd);
    
    sumCorrROI = matR_SU_all_moviemask'*matClusterInd;
    meanCorrROI = sumCorrROI./repmat(nVoxPerROI, size(sumCorrROI, 1), 1); 
    [sortMeanCorrROI, indCell] = sortrows(meanCorrROI, -1.*(clusterOrder_arbitrary{iK}(1)));

    
    figCorrMat_voxCluster_sorted = figure;
    set(figCorrMat_voxCluster_sorted, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 300 1200])
    imagesc(sortMeanCorrROI(:, clusterOrder_arbitrary{iK}))
    set(gca, 'CLim', [-1 1].*0.6)
    set(gca, 'XTick', 1:targetK, 'XTickLabel', clusterOrder_arbitrary{iK}, 'YTick', [])
    colorbar('NorthOutside')
    title(sprintf('Correlation between neurons and voxels: %d ROIs', targetK))
    xlabel('Voxel clusters')
    ylabel('Cells')
    print(figCorrMat_voxCluster_sorted, fullfile(dirFig, sprintf('CorrMatrix_%s%s_%dROI_sorted', cell2mat(setNameSubjNeural), nameSubjBOLD, targetK)),...
        '-dtiff', '-r200')
end






