% analClustering_multiplePatches.m
% 
% 2019/06/18 Check the clustering results: based on whole-brain
% correlations


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
addpath(fullfile(dirProjects, 'parksh/_toolbox/afni_matlab'))
addpath(fullfile(dirProjects, 'parksh/_toolbox/hslcolormap'))

% Set directories
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was'};
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh/_macaque');
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

%% Load data
% 1) Clustering based on corr maps
load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)));
% load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_pcares_masked.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)));  
% % load(fullfile(dirDataBOLD, sprintf('Clustering_%s%sMovie123_new_masked_probability_critCorr1.mat', cell2mat(setNameSubjNeural), nameSubjBOLD)))  

% load the masks (movie-driven & brain-only mask)
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');


%% Load the corr map and select valid channel: subject by subject
numSubject = size(setNameSubjNeural, 2);
matR_SU_all_moviemask = [];

for iSubj = 1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');
%     load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_pcares_masked.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');

    matR_SU_valid = matR_SU(:, paramClustering(iSubj).validChanIndex);
    matR_SU_all_moviemask = cat(2, matR_SU_all_moviemask, matR_SU_valid);
    
     clear matR_SU matR_SU_valid
end

% [nx ny nz] = size(movieDrivenAmp.mask_amp1);
% nVox = nx*ny*nz;

% Select voxels: voxels that showed correlation higher than 0.3 with any one
% of the neurons
critCorr = 0.3;
critNumNeuron = ceil(size(matR_SU_all_moviemask, 2).*0.05); % 5% of the population %6; %13;
matValidVox = abs(matR_SU_all_moviemask)>critCorr;
locValidVox = sum(matValidVox, 2)>critNumNeuron;
matR_SU_all_moviemask_valid  = matR_SU_all_moviemask(locValidVox, :);


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
    
    matWSS(:,iK) = sum(Clustering_moviemask_valid.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
%     matExpVar(iK,1) = sum(tExpVar);

end

totalSS = Clustering_moviemask_valid.totalSS_SU;
% [a, c, totalSS] = kmeans(matR_SU_moviemask', 1); %kmeans(matR_SU', 1);
betweenSS = totalSS-matWSS;
% totalVar = totalSS/(2*size(matR_SU,2)); %totalD/(2*size(matR_SU,2));

propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;


figure;
plot(setK, propExplained', 'ko-'); hold on
xlabel('Number of cluster (K)')
ylabel('Explained variance (%)')
set(gca, 'XTick', setK)
% set(gca, 'YTick', 0.55:0.05:0.9, 'YTicklabel', 55:5:90)
% curK=5;
% plot(setK, propExplained', 'ko-', 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 8); hold on
% plot(curK, propExplained(curK-1), 'ko', 'LineWidth', 2, 'MarkerFaceColor', 'k', 'MarkerSize', 8)
% xlim([2 12])
% ylim([0.55 0.9])
% set(gca, 'TickDir', 'out', 'LineWidth', 2, 'Box', 'off', 'TickLength', [.025 .05])
% set(gca, 'YTick', 0.55:.1:.9)

% plot Clustering results
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iK=1:length(setK)
    plot(iK+1, Clustering_moviemask_valid.resultKMeans(iK).SU_sumD, 'ko', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
end
xlim([1 setK(end)])
set(gca, 'LineWidth', 2, 'FontSize', 12)
box off
title('Clustering of single units')
xlabel('Number of cluster')
ylabel('Within-cluster distance')

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
for iK=1:length(setK)
    plot(iK+1, mean(Clustering_moviemask_valid.resultKMeans(iK).SU_totalD), 'ko', 'MarkerSize', 10, 'LineWidth', 2)
    hold on
end
xlim([1 setK(end)])
set(gca, 'LineWidth', 2, 'FontSize', 12)
box off
title('Clustering of single units')
xlabel('Number of cluster')
ylabel('Within-cluster distance')

% % save
% print(fig3b, fullfile(dirFig, 'figRev3B_VarExplained'), '-depsc')

%% MDS plot showing K-means clustering results
D = pdist(matR_SU_all_moviemask', 'euclidean');
[Y2,stress,disparities] = mdscale(D,2);
[Y3,stress,disparities] = mdscale(D,3);

% curK=7;
% locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
% locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
% [sortedClust, indSortChan] = sort(Clustering_moviemask_valid.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));

matIndClust_SU = cat(2, Clustering_moviemask_valid.resultKMeans.SU_indCluster); %cat(2, Clustering.resultKMeans.SU_indCluster);
curK = 14; % 6; %12; %6; %4; %7; %
[sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));
indNewCluster = 1:curK; %[4 2 3 6 7 5 1]; % [4 2 5 1 7 6 3]; %[1 2 3 4 5]; %[1 3 6 4 2 5]; %[2 4 1 3]; %[1 2 3 4 5]; %[3 4 2 5 1]; %[1 2 3 4 5 6 7]; %[4 1 6 3 5 2 7]; % cluster #4 is now cluster 1

% % index for different face patches: AF = 1; AM = 2; peri-AM = 3;
% for iSubj = 1:length(paramClustering)
%     if iSubj<5 % AF
%         catIndArea = cat(1, catIndArea, ones(size(paramClustering(iSubj).validChanIndex).*1)); 
%     elseif ismember(iSubj, [5 8]) % Matcha and Wasabi
%         catIndArea = cat(1, catIndArea, ones(size(paramClustering(iSubj).validChanIndex).*2)); 
%     elseif iSubj == 6 % Dango
%         catIndArea = cat(1, catIndArea, ones(size(paramClustering(iSubj).validChanIndex).*3));
%     elseif iSubj == 7 % Mochi, has AF and AM
%         aa= regexp(cellstr(paramClustering(iSubj).validChanID), '[A].', 'match');
        

% indMonkey = [];
for iC = 1:length(paramClustering)
    curTS(iC).validChanID = strcat(cellstr(paramClustering(iC).validChanID), paramClustering(iC).nameSubj); %cellstr(paramClustering(iC).validChanID);
%     indMonkey = cat(1, indMonkey, ones(size(paramClustering(iC).validChanIndex)).*iC);
end
catChanID = cat(1, curTS.validChanID);
catSubjID = cat(1, paramClustering.validChan_subjID);
catAreaID = floor(catSubjID./10);
% catChanID = cat(1, paramClustering.validChanID);
% indMonkey = str2num(catChanID(:,5));

% Directory for saving figures as graphic files
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');


% Plotting parameters
cMap = cool(3);
% % cMap = [0 0 0; 230 159 0; 86 180 233; 0 158 115; 240 228 66; 0 114 178; 213 94 0; 204 121 167]./255;
% cMap = [228	26	28;
% 55	126	184;
% 77	175	74;
% 152	78	163;
% 255	127	0;
% 255 217 47; % dark yellow %255	255	255; %white %255	255	51; % yellow was too similar to another yellow in mat2
% 166	86	40;
% 247	129	191;
% ]./255;
marker = {'o', '^', 'v', '<', '>', 'square', 'diamond', '*'}; %{'o', '^', 'square', 'diamond'}; % {'o','*', 'x', 's', 'd', '+', '^'};
% orderIndNewClust = [3 2 7 6 5 1 4];
% cMap_newclust = cMap(orderIndNewClust,:);

%
fig3a2=figure;
set(fig3a2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 600 700 700])
for iC=1:curK
    
    iK = indNewCluster(iC);
    
    curIndCells = indSortChan(sortedClust==iK); %find(sortedClust==iK); %resultProbClustering(iK).validIndCells;
    combCells = nchoosek(curIndCells,2);

    figure(fig3a2);

%     line([Y2(combCells(:,1), 1) Y2(combCells(:,2),1)]', [Y2(combCells(:,1), 2) Y2(combCells(:,2), 2)]', 'Color', cMap(iC,:), 'LineWidth', 2);
%     hold on;
    for iCell = 1:length(curIndCells)
        plot(Y2(curIndCells(iCell), 1), Y2(curIndCells(iCell), 2), 'o-',... % 'Marker', marker{indMonkey(curIndCells(iCell))},...
            'MarkerSize', 8, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(catAreaID(curIndCells(iCell)),:), 'LineWidth', 2); %cMap(iC,:), 'LineWidth', 2)
        text(Y2(curIndCells(iCell),1)-.5, Y2(curIndCells(iCell),2)+.2, catChanID(curIndCells(iCell),:))
        
%         plot3(Y3(curIndCells(iCell),1), Y3(curIndCells(iCell),2), Y3(curIndCells(iCell),3), 'o-',... %'Marker', marker{indMonkey(curIndCells(iCell))}, ...
%             'LineWidth', 2, 'MarkerSize', 8,...
%             'MarkerEdgeColor', 'k', 'MarkerFaceColor', cMap(catAreaID(curIndCells(iCell)),:)); % cMap(iC,:));
%         text(Y3(curIndCells(iCell),1), Y3(curIndCells(iCell),2), Y3(curIndCells(iCell),3), catChanID(curIndCells(iCell),:))
        hold on
    end    
%     figure(fig3a2); 
% %     curChan = indSortChan(sortedClust==iK);
%     p = plot(Y2(curIndCells, 1), Y2(curIndCells, 2), 'o','LineWidth', 2, 'MarkerSize', 10,...
%         'MarkerEdgeColor','k', 'MarkerFaceColor',cMap(iC,:), 'Color', cMap(iC,:)); %,...
%         %'Marker', marker{indMonkey(resultProbClustering(iC).validIndCells)});
%     hold on;
% %     plot(Y2(resultProbClustering(iC).validIndCells,1), Y2(resultProbClustering(iC).validIndCells,2), 'o','LineWidth', 2, 'MarkerSize', 10,...
% %         'MarkerEdgeColor','k', 'MarkerFaceColor', cMap(iC,:), 'Color', cMap(iC,:));
% %     text(Y2(curChan,1)+1, Y2(curChan,2), catChanID(curChan,:)) %paramCorr.validChanID(curChan,:))
%     hold on;
%     xlim([-6 8])
%     ylim([-8 8])
%     zlim([-4 4])
input ('')
end

axis square
xlim([-32 35]) %xlim([-35 35]) %xlim([-25 30])
ylim([-10 11]) %ylim([-12 12])
set(gca, 'XTick', [], 'YTick', [], 'LineWidth', 2,  'Box', 'off')


%% Check the clustering
% Visualization of similarity of maps within each cluster 
for iC = 1:length(paramClustering)
    temp(iC).validChanID = strcat(cellstr(paramClustering(iC).validChanID), paramClustering(iC).nameSubj); %cellstr(paramClustering(iC).validChanID);
%     indMonkey = cat(1, indMonkey, ones(size(paramClustering(iC).validChanIndex)).*iC);
end
catChanID = cat(1, temp.validChanID);
catSubjID = cat(1, paramClustering.validChan_subjID);
catAreaID = floor(catSubjID./10);

matIndClust_SU = cat(2, Clustering_moviemask_valid.resultKMeans.SU_indCluster); %cat(2, Clustering.resultKMeans.SU_indCluster);
curK = 14; %7; % 14; % 6; %12; %6; %4; %7; %
[sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));

figure
imagesc(matR_SU_all_moviemask)
set(gca, 'CLim', [-1 1].*0.5)
figure
imagesc(matR_SU_all_moviemask(:, indSortChan))
locClust = find(diff(sortedClust)>0);
set(gca, 'XTick', find(diff(sortedClust)>0))

% Visualization of area composition for each cluster
figure;
imagesc(catAreaID(indSortChan)); % Blue: AF, Green: AM, Red: AM+
set(gca, 'YTick', find(diff(sortedClust)>0)) 


%% Plot the mean time series for each cluster
load('/procdata/parksh/_macaque/matSDF_Movie123_allCells.mat')

matFR_TR = cat(2, matSDF.matFR_TR);
matNeuralRGR = cat(2, matSDF.matNeuralRGR);
matFR_SU = cat(2, matSDF.matFR_SU_10hz);
matFR_SU_norm = cat(2, matSDF.matFR_SU_10hz_norm);

% matTS = NaN(size(matFR_TR, 1), 3

for iC = 1:length(paramClustering)
    temp(iC).validChanID = strcat(cellstr(paramClustering(iC).validChanID), paramClustering(iC).nameSubj); %cellstr(paramClustering(iC).validChanID);
%     indMonkey = cat(1, indMonkey, ones(size(paramClustering(iC).validChanIndex)).*iC);
end
catChanID = cat(1, temp.validChanID);
catSubjID = cat(1, paramClustering.validChan_subjID);
catAreaID = floor(catSubjID./10);


matIndClust_SU = cat(2, Clustering_moviemask_valid.resultKMeans.SU_indCluster); % cat(2, ClusteringSDF_TR.resultKMeans.SU_indCluster); %%cat(2, Clustering.resultKMeans.SU_indCluster);
curK = 10; %14; %7; %14; %7; % 14; % 6; %12; %6; %4; %7; %
[sortedClust, indSortChan]=sort(matIndClust_SU(:,curK-1));

figTS = figure;
set(figTS,  'Color', 'w', 'PaperPositionMode', 'auto')
cMap = eye(3).*0.5;
matTS=[];
for iK = 1:curK
    curIndCells = indSortChan(sortedClust==iK);    
    curArea = catAreaID(curIndCells);
    
    curTS = zscore(matFR_TR(:,curIndCells));
%     curTS = matNeuralRGR(:,curIndCells);
%     curTS = matFR_SU_norm(:,curIndCells);
    
    meanTS = []; steTS = []; 
    for iA = 1:3
        meanTS(:, iA) = nanmean(curTS(:, curArea==iA), 2);
        steTS(:, iA) = nanstd(curTS(:, curArea==iA), 0, 2)./sum(curArea==iA);
    end
%     tempCat = cat(3, meanTS+steTS, meanTS-steTS);
%     bb = reshape(aa, 375*3, 2);
    figure(figTS); clf;
%     subplot(ceil(curK/2), 2, iK);
    for iA = 1:3
        figure(figTS); 
        line(repmat(1:size(meanTS, 1), 2, 1), cat(2, meanTS(:,iA)+steTS(:,iA), meanTS(:,iA)-steTS(:,iA))', 'Color', cMap(iA, :), 'LineWidth', 2);
        hold on;
        PL(iA) = plot(1:size(meanTS, 1), meanTS(:, iA), 'o-', 'Color', cMap(iA, :), 'LineWidth', 2, ...
             'MarkerEdgeColor', cMap(iA, :), 'MarkerFaceColor', 'w');
%          PL(iA) = plot(1:size(meanTS, 1), meanTS(:, iA), '-', 'Color', cMap(iA, :), 'LineWidth', 2);
        hold on;        
    end
    
    legend(PL, sprintf('AF: n = %d', sum(curArea==1)), sprintf('AM: n = %d', sum(curArea==2)), sprintf('AM+: n = %d', sum(curArea==3)))
    title(sprintf('Averaged time series of neurons from each area in each cluster: Cluster # %d / %d', iK, curK))
    input('')

    matTS = cat(2, matTS, meanTS);
end
    
        
%% Compare different clustering
% Load and retrieve necessary struct
load('/procdata/parksh/_macaque/Art/Clustering_TorRhoSigSpiMatDanMocWasArtMovie123_pcares_masked.mat');
Clustering_moviemask_pcares = Clustering_moviemask;
clear Clustering_moviemask;
load('/procdata/parksh/_macaque/Art/Clustering_TorRhoSigSpiMatDanMocWasArtMovie123_new_masked_critCorr1.mat')
Clustering_moviemask_valid_orgCorrMap = Clustering_moviemask_valid;
clear Clustering_moviemask_valid;

% Sort out cells based on particular K-means clustering 
matIndClust_SU = cat(2, Clustering_moviemask_valid_orgCorrMap.resultKMeans.SU_indCluster); % based on corr map
matIndClust_SU_pcares = cat(2, Clustering_moviemask_pcares.resultKMeans.SU_indCluster); % based on pca-res corr map

sortTargetK = 7;
[sortedClust, indSortChan]=sort(matIndClust_SU(:,sortTargetK-1));
[sortedClustSDF, indSortChanSDF]=sort(matIndClust_SU_pcares(:,sortTargetK-1));


% Comparison figure
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
image(cat(2, matIndClust_SU(indSortChan,sortTargetK-1),... % matIndClustSDF_SU(indSortChan, sortTargetK-1), ...
    matIndClustSDFnorm_SU(indSortChan, sortTargetK-1), matIndClustSDFTR_SU(indSortChan, sortTargetK-1), ...
    matIndClustSDFMION_SU(indSortChan, sortTargetK-1), matIndClustSDFRGR_SU(indSortChan, sortTargetK-1)))
colormap(lines)
% h=colorbar;
% set(h, 'YLim', [0 sortTargetK]+0.5)
% set(h, 'YTick', [])
set(gca, 'XTick', [])
set(gca, 'YTick', 1:48, 'YTickLabel', paramCorr.validChanID(indSortChan,:))
set(gca, 'FontSize', 12)
title('Results of K-means clustering of cells')
set(gca, 'XTick', 1:5, 'XTickLabel', {'Corr maps', 'Time Series (zscore)', 'Time Series (TR)', 'Time Series (MION)', 'Movie feature correlation'}); % {'Corr maps', 'Time Series (fine)', 'Time Series (zscore)', 'Time Series (TR)', 'Time Series (MION)'})
xlabel('Basis of clustering')
ylabel('Cells')
    


       



% dirDataBOLD = '/procdata/parksh/_macaque/Art';
% dirROI = fullfile(dirDataBOLD, 'ROIs');
% % if ~exist(dirROI)
% %     visROIs=[]; faceROIs=[];
% % else
% load(fullfile(dirROI, 'e66_faceROIs_L.mat'))
% load(fullfile(dirROI, 'e66_faceROIs2.mat'))
% 
% faceROIs = struct([]);
% faceROIs = cat(2, e66_faceROIs2, e66_faceROIs_L);
% 
% clear e66*
% 
% % Get ROI names & coordinates in EPI coords 
% resizeFactor = [3 3 3]; % DSP.proc.params3d.res./ROIs(1,1).params.res; % calculate the scaling factor
% 
% for iR = 1:length(faceROIs) 
% nameROI = faceROIs(iR).name; %(strfind(faceROIs(iR).name, '_')+1:end); %faceROIs(iFR).name(strfind(faceROIs(iFR).name, '_')+1:end);
% voxROI=decimate3D(faceROIs(iR).vol3D, resizeFactor, .25); %decimate3D(faceROIs(iFR).vol3D, resizeFactor, .25); % turn anat_res ROIs into func_res ROIs
% [a, b, c] = ind2sub(size(voxROI), find(voxROI==1)); % Get indices of AF voxels in EPI 3D coords
% indVox_ROI = [a b c];
% indVox_ROI_sub = sub2ind(size(voxROI), a, b, c);

