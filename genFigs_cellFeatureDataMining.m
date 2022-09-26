% genFigs_cellFeatureDataMining.m
%
% 2022/09/05 SHP
% Resurrection of the old attempt
% 2017/10/09 SHP
% Wondering around to find out something in cells' responses to movies
% in relation to the features or eye gaze (see also
% genFigs_eyeGazeDataMining.m)

clear all;

ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/NIFVAULT/projects';
    dirProcdata = '/Volumes/NIFVAULT/procdata';
    dirDataHome = fullfile(dirProcdata, 'parksh/_macaque');
    dirLibrary = '/Volumes/NIFVAULT/library';
    %         addpath(fullfile(dirLibrary, 'matlab_utils'));
%     dirDataNeural = fullfile(dirDataHome, 'Spi');
else % on virtual machine
    dirProjects = '/nifvault/projects';
    dirProcdata = '/nifvault/procdata';
    dirDataHome = fullfile(dirProcdata, 'parksh/_macaque');
    dirLibrary = '/nifvault/library';
    %         addpath(fullfile(dirLibrary, 'matlab_utils'));
%     dirDataNeural = fullfile(dirDataHome, 'Spi');
end

%% Load the cell time series & movie regressors
load(fullfile(dirDataHome, 'matSDF_Movie123_30fps.mat'))
load(fullfile(dirDataHome, 'movieRegressors_4fps_selected.mat')); %'fullMovieRegressors_30fps.mat'))

tempS = cat(1, matSDF.FR_30fps);
indFP = ~contains(cat(1, matSDF.validChanID), 'NFP');
FR_30fps_FP = tempS(indFP, :);
FR_30fps_NFP = tempS(~indFP, :);

matTS = []; matRGR = [];
for iMovie = 1:3 %iMovie = 3;    
    matTS = cat(1, matTS, cat(2, FR_30fps_FP(:,iMovie).mnFR));
    matRGR = cat(1, matRGR, resample(fullRGR4fps(iMovie).regressors, 30, 4));
end

% Sort the cell TS along the area of the cells
catSubjID = cat(1, matSDF.validChan_subjID);
catSubjID_FP = catSubjID(indFP);
[a, b] = sort(catSubjID_FP);
matTS = matTS(:, b);
catSubjID_FP = catSubjID_FP(b);
catAreaID_FP = floor(catSubjID_FP./10);

%% Correlation between individual neurons and movie regressors
R_SUmovieRGRfull=NaN(size(matRGR,2), size(matTS,2));
for iRGR = 1:size(matRGR,2)
    r_c=[]; r_su=[];
    
%     % averaged TS in each cluster
%     r_c = corr(matRGRfull(:,iRGR), meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');
    % single unit TS
    r_su = corr(matRGR(:,iRGR), matTS, 'rows', 'complete', 'type', 'Spearman');
    
%     R_ClusterMovieRGRfull(iRGR, :) = r_c;
    R_SUmovieRGRfull(iRGR, :) = r_su;    
    
end

% subset of regressors
indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
R_SUmovieRGRvalid = R_SUmovieRGRfull(indValidRGR,:);
varnamesvalid = fullRGR4fps(1).features(indValidRGR);
% R_ClusterMovieRGRvalid = R_ClusterMovieRGRfull(indValidRGR,:);



%% correlation between principal components and movie regressors
%% PCA of face-selective neurons
load('/nifvault/procdata/parksh/_macaque/multipleFP_fsi.mat')
locFaceCell =  find(fsi.matFSI(:,2)>0); % find(abs(fsi.matFSI(:,1))>0.33);
locNonFaceCell = find(fsi.matFSI(:,2)<0);

matTS_norm = zscore(matTS); %zscore(matTS(:, locFaceCell));
[coeff, score, latent, tsq, explained] = pca(matTS_norm');
[sortedScore, indSortChan] = sort(score(:,1));

% Correlation 
R_PCmovieRGRfull=NaN(size(matRGR,2), 10);
for iRGR = 1:size(matRGR,2)
    r_c=[]; r_pc=[];
    
%     % averaged TS in each cluster
%     r_c = corr(matRGRfull(:,iRGR), meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');
    % single unit TS
    r_pc = corr(matRGR(:,iRGR), coeff(:, 1:10), 'rows', 'complete', 'type', 'Spearman');
    
%     R_ClusterMovieRGRfull(iRGR, :) = r_c;
    R_PCmovieRGRfull(iRGR, :) = r_pc;    
    
end

varnamesfull = fullRGR4fps(1).features;
% reorder the full regressors 
indReorderRGR = [1 2 9 10 6 7 8 3 11:18 4 5 21 27 28 20 29 30 22 31 32:36 23 24 25 26 19];
% catChanID = cat(1, paramClustering.validChanID);

% 1. Correlation matrix between SUs and regressors
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1600 600])
imagesc(R_PCmovieRGRfull(indReorderRGR, :)) % indSortChan_new));


%% Plot
varnamesfull = fullRGR4fps(1).features;
% reorder the full regressors 
indReorderRGR = [1 2 9 10 6 7 8 3 11:18 4 5 21 27 28 20 29 30 22 31 32:36 23 24 25 26 19];
% catChanID = cat(1, paramClustering.validChanID);

% 1. Correlation matrix between SUs and regressors
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1600 600])
imagesc(R_SUmovieRGRfull(indReorderRGR, :)) % indSortChan_new));
xPos = [find(diff(catAreaID_FP)>0)+0.5; size(R_SUmovieRGRfull, 2)];
line([xPos xPos]', repmat([0; 37], 1, length(xPos)), 'Color', 'k') 
cval = 0.5;
set(gca, 'CLim', [-1 1].*cval, 'FontSize', 12)
set(gca, 'XTick', xPos, 'XTickLabel', [])
set(gca, 'YTick', 1:length(varnamesfull), 'YTickLabel', varnamesfull(indReorderRGR,:))
% make blue-white-red colorbar

cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end
endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));

figure(gcf)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;
title('Correlation between SU time series and feature time series')

% print(gcf, fullfile(dirFig, 'Corr_SUMovieRGR_movie123_ClusteringCorrMapTorRhoSigSpiArt'), '-dtiff', '-r150')


%% fewer selected regressors
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1600 600])
imagesc(R_SUmovieRGRvalid) % indSortChan_new));
xPos = [find(diff(catAreaID_FP)>0)+0.5; size(R_SUmovieRGRvalid, 2)];
line([xPos xPos]', repmat([0; length(varnamesvalid)]+1, 1, length(xPos)), 'Color', 'k') 
cval = 0.5;
set(gca, 'CLim', [-1 1].*cval, 'FontSize', 12)
set(gca, 'XTick', xPos, 'XTickLabel', [])
set(gca, 'YTick', 1:length(varnamesvalid), 'YTickLabel', varnamesvalid)
% make blue-white-red colorbar

cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end
endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));

figure(gcf)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;
title('Correlation between SU time series and feature time series')

%% face vs. non-face cells?
load('/nifvault/procdata/parksh/_macaque/multipleFP_fsi.mat')
locFaceCell =  find(fsi.matFSI(:,2)>0); % find(abs(fsi.matFSI(:,1))>0.33);
locNonFaceCell = find(fsi.matFSI(:,2)<0);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1600 600])
imagesc(R_SUmovieRGRvalid(:, [locFaceCell; locNonFaceCell])) % indSortChan_new));
xPos = [length(locFaceCell)+0.5; length(locFaceCell) + length(locNonFaceCell)];
line([xPos xPos]', repmat([0; length(varnamesvalid)+1], 1, length(xPos)), 'Color', 'k') 
cval = 0.5;
set(gca, 'CLim', [-1 1].*cval, 'FontSize', 12)
set(gca, 'XTick', xPos, 'XTickLabel', [])
set(gca, 'YTick', 1:length(varnamesvalid), 'YTickLabel', varnamesvalid)
% make blue-white-red colorbar

cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end
endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));

figure(gcf)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;
title('Correlation between SU time series and feature time series')



%% order the face cells according to fMRI-map based clustering
load('/nifvault/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_faceselective_Movie123_probability.mat', 'Clustering_brainmask', 'param*') 
% load('/nifvault/procdata/parksh/_macaque/Art/Clustering_CorrMap_4FPs_Movie123_probability.mat', 'Clustering_brainmask', 'param*') 
% load(sprintf('/procdata/parksh/_macaque/CorrMap_SU_AllCells%s_corticalFPMerged.mat', nameSubjBOLD), 'corrMap_merged_FP'); %, 'info*', 'corrMap_Area', 'corrMap_merged');

setK = 2:20; %paramClustering_global.setK; %Clustering.setK;

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS(:,iK) = sum(Clustering_brainmask.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
end
totalSS = Clustering_brainmask.totalSS_SU;
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;

% Get the clustering results from voxel-based clustering
curK = 6; %10; %13; % 9; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

clusterS = struct([]);
% ordering of cells: ML-AF-AM-AAM, while maintaining the grouping
for iK = 1:curK
    clusterS(iK).indSortChan_org = indSortChan(sortedClust==iK);
    clusterS(iK).sortedClust_org = sortedClust(sortedClust==iK);
    clusterS(iK).numCells = sum(sortedClust==iK);
    locML = find(Clustering_brainmask.infoCells.catAreaID(clusterS(iK).indSortChan_org)==4);
    if isempty(locML)
        clusterS(iK).indSortChan_reorder_MLfirst = clusterS(iK).indSortChan_org;
    else
        clusterS(iK).indSortChan_reorder_MLfirst = cat(1, clusterS(iK).indSortChan_org(locML), clusterS(iK).indSortChan_org(1:locML(1)-1));
    end
end
[~,reorderCluster] = sort(cat(1, clusterS.numCells), 'descend');
% reorderCluster = [8 1 9 5 3 7 4 10 6 2]; % ORder cell groups based on the number of neurons in each group (from largest to smallest)
indSortChan_reorder = cat(1, clusterS(reorderCluster).indSortChan_org);
sortedClust_reorder = cat(1, clusterS(reorderCluster).sortedClust_org);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1600 600])
imagesc(R_SUmovieRGRvalid(:, locFaceCell(indSortChan_reorder))) % indSortChan_new));
xPos = [find(abs(diff(sortedClust_reorder))>0); length(locFaceCell)];
line([xPos xPos]', repmat([0; length(varnamesvalid)+1], 1, length(xPos)), 'Color', 'k') 
cval = 0.5;
set(gca, 'CLim', [-1 1].*cval, 'FontSize', 12)
set(gca, 'XTick', xPos, 'XTickLabel', [])
set(gca, 'YTick', 1:length(varnamesvalid), 'YTickLabel', varnamesvalid)
% make blue-white-red colorbar

cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end
endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));

figure(gcf)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;
title('Correlation between SU time series and feature time series')


% % 1. Correlation matrix between clusters and regressors
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 650 880])
% imagesc(R_ClusterMovieRGRfull(indReorderRGR, :));
% colorbar;
% set(gca, 'CLim', [-.6 .6], 'FontSize', 12)
% set(gca, 'YTick', 1:length(varnamesfull), 'YTickLabel', varnamesfull(indReorderRGR,:))
% xlabel('Cell Group')
% title('Correlation between cluster time series and feature time series')
% % make blue-white-red colorbar
% cval = 0.6;
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
% set(gca, 'TickDir', 'out')
% box off
% c=colorbar;
% 
% print(gcf, fullfile(dirFig, 'Corr_ClusterMovieRGR_movie123_ClusteringCorrMapTorRhoSigSpiArt'), '-dtiff', '-r150')


% indValidRGR = [1     2     3     4     5     6     7     8     9    10    11    12    13    14    15    16    17    18 ...
%     19    20    21    22    23    25    26    27    28    29    30    31    32    33    34    36    38    39 ...
%     40    41    42    46    47    52    55    57    58    65    66    67    68    69    70    71]';
% varnamesvalid = fullRGR30fps(1).features(indValidRGR);
% 
% iMovie = 3;
% 
% catTS = cat(2, FR_30fps_FP(:,iMovie).mnFR);
% catValidRGR = fullRGR30fps(iMovie).regressors(:, indValidRGR);
% 
% R_SUmovieRGRfull = NaN(length(indValidRGR), size(catTS, 2));
% for iRGR = 1:length(indValidRGR) %size(fullRGR30fps(iMovie).regressors, 2)
%     r_su = [];
%     r_su = corr(catValidRGR(:, iRGR), catTS, 'rows', 'complete', 'type', 'Spearman');
%     
%     R_SUmovieRGRfull(iRGR, :) = r_su;
% end
% 
% figure
% imagesc(R_SUmovieRGRfull);
% colormap(jet)
% set(gca, 'CLim', [-1 1].*0.5)
% set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', varnamesvalid)
% colorbar


% setMovie = [1 2 3];
% load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr') % we only need ID of valid channels
% 
% % Prepare neural responses
% % First, get the cell response in 30fps time resolution 
% FR_dT30 = createCellRegressor_indMov_discreteTime(dirDataNeural, cellstr(paramCorr.validChanID),...
%     setMovie, 1/30); 
% % concatenate across movies
% matFR=[];
% for iUnit = 1:size(FR_dT30,1)
%     tempFR = cat(1, FR_dT30(iUnit, :).mnFR);
%     matFR(:,iUnit) = tempFR;
% end
% 
% % subset of regressors
% indValidRGR = [1 2 6 7 3 21 20 32 22 31 25]; %[1, 3, 9, 20, 21, 22, 25]; 
% % 1: 'Luminance', 2: 'Contrast', 6: Low spatial Frequency 7: High spatial frequencty 3: 'Motion (Speed)', 
% % 21: 'One face', 20: 'Number of faces', 32: 'Face size', 22: 'Body parts', 31: 'Hands', 25: 'Any animal'
% % matRGRvalid = matRGRfull(:,indValidRGR);
% varnamesvalid = varnamesfull(indValidRGR);
% 
% % Compute correlation between neural TS and feature TS
% R_ClusterMovieRGRfull=NaN(size(matRGRfull,2), size(meanFRCluster4fps,2));
% R_SUmovieRGRfull=NaN(size(matRGRfull,2), size(matFR4fps,2));
% for iRGR = 1:size(matRGRfull,2)
%     r_c=[]; r_su=[];
%     
%     % averaged TS in each cluster
%     r_c = corr(matRGRfull(:,iRGR), meanFRCluster4fps, 'rows', 'complete', 'type', 'Spearman');
%     % single unit TS
%     r_su = corr(matRGRfull(:,iRGR), matFR4fps, 'rows', 'complete', 'type', 'Spearman');
%     
%     R_ClusterMovieRGRfull(iRGR, :) = r_c;
%     R_SUmovieRGRfull(iRGR, :) = r_su;    
%     
% end
% 
% R_SUmovieRGRvalid = R_SUmovieRGRfull(indValidRGR,:);
% R_ClusterMovieRGRvalid = R_ClusterMovieRGRfull(indValidRGR,:);


% % setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
% setMovie = [1 2 3];
% % nt = 375;
% % 
% % % MION function
% % TR=2.4;
% % k = gampdf([-40:TR:40],4,2);
% 
% matSDF = struct([]);
% curUnit = 0;
% for iSubj = 1:length(setNameSubjNeural)
%     
%     % Load the parameter
%     nameSubjNeural = setNameSubjNeural{iSubj};
%     dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
%     filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
%     load(fullfile(dirDataNeural, filenameNeural), 'paramSDF')
%     fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
%         
%     switch lower(nameSubjNeural)
%         case 'spi'
%             excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
%             validC = setdiff(1:length(paramSDF.setCellIDs), excChanIndex)';
%         otherwise
%             [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
%             validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
%     end   
%     matSDF(iSubj).setCellIDs = paramSDF.setCellIDs(validC);
%     
%     indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
%     % 3. Fine temporal resolution
%     clear FR_dTfine
%     FR_dTfine = createCellRegressor_indMov_discreteTime(dirDataNeural, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
%         setMovie, 0.1); % in 10Hz (number of spikes for every 100ms)
%     % concatenate across movies
% %     matFR_SU=cell(135, );
%     for iUnit = 1:size(FR_dTfine,1)
%         curUnit = curUnit + 1;
%         tempFR = cat(2, FR_dTfine(iUnit, :).matFR); 
%         matFR_SU(curUnit, :) = tempFR; 
%         
%         %cat(1, FR_dTfine(iUnit, :).mnFR);
%         %matFR_SU(:,iUnit) = tempFR;
%     end
% %     matSDF(iSubj).matFR_SU = matFR_SU;
%     
% end
% 
% % for iCell = 1:size(matFR_SU, 1)
% % for iMov = 1:size(matFR_SU, 2)
% % %     numTrial(iCell, iMov) = size(FR_dTfine(iCell, iMov).matFR{1}, 2);
% % avgSTD(iCell, iMov) = mean(std(matFR_SU{iCell, iMov}'));
% % end
% % end
% 
% for iCell = 1:size(matFR_SU, 1)
% for iMov = 1:size(matFR_SU, 2)
% tempR = corrcoef(matFR_SU{iCell, iMov});
% avgR(iCell, iMov) = mean(tempR(tril(tempR,-1)>0));
% end
% end
% 
% % 1) Clustering based on corr maps
% load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  
% 
% 
% %% 
% setK = paramClustering_global.setK; %Clustering.setK;
% 
% matWSS_corrMap=[];
% matExpVar=[];
% for iK = 1:length(setK)
%     curK = setK(iK);
%     matWSS_corrMap(:,iK) = sum(Clustering_moviemask_valid.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
% end
% 
% totalSS_corrMap = Clustering_moviemask_valid.totalSS_SU;
% % betweenSS_corrMap = totalSS_corrMap-matWSS_corrMap;
% propExplained_corrMap = (totalSS_corrMap-matWSS_corrMap)./totalSS_corrMap; %matExpVar./totalSS;
% 
% 
% %% Compare different clustering
% K=7;
% locMode_corrMap = find(propExplained_corrMap(:,K-1)==mode(propExplained_corrMap(:,K-1)));
% indClust_SU = Clustering_moviemask_valid.resultKMeans(K-1).SU_indCluster(:, locMode_corrMap(1)); % based on corr map
% [sortedClust, indSortChan]=sort(indClust_SU);
% 
% oldIndCluster = [1 2 5 6 3 7 4]; % [3 4 5 2 1]; % [4 1 6 3 5 2 7]; 
% indSortChan_new = [];
% for iC = 1:K %7
%     curC = oldIndCluster(iC);
%     tempind = indSortChan(sortedClust==curC);
%     indSortChan_new = cat(1, indSortChan_new, tempind);
% end
% sortedClust_new = indClust_SU(indSortChan_new);
% 
% xPos = find(abs(diff(sortedClust_new))>0)+0.5;
% line([xPos xPos]', repmat([0; 37], 1, 6), 'Color', 'k') 
%  
