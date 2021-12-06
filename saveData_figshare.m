% saveData_figshare.m
%
% 2021/12/06 SHP
% Generate the dataset for sharing via figshare
% in relation to Park et al., 2021

clear all;

%% Settings
flagBiowulf = 1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        directory.projects = '/Volumes/PROJECTS';
        directory.procdata = '/Volumes/PROCDATA';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/Volumes/LIBRARY';
        addpath(fullfile(directory.library, 'matlab_utils'));
    else % on virtual machine
        directory.projects = '/nifvault/NIFVAULT/projects';
        directory.procdata = '/procdata';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/library';
        addpath(fullfile(directory.library, 'matlab_utils'));
    end
end

%%
load(fullfile(directory.dataHome, 'Art/Clustering_CorrMap_4FPs_Movie123_ArtRHROI_set01_probability.mat'), 'Clustering_meanROI')
load(fullfile(directory.dataHome, 'Art/Clustering_CorrMap_4FPs_faceselective_Movie123_probability.mat'), 'Clustering_brainmask')
load(fullfile(directory.dataHome, 'multipleFP_fsi.mat'))

%%
locFaceCell =  find(fsi.matFSI(:,1)>0.33); 
matR = Clustering_meanROI.matR(locFaceCell, :);
catAreaID = Clustering_meanROI.catAreaID(locFaceCell);
catChanID = Clustering_meanROI.catChanID(locFaceCell);
setArea = Clustering_meanROI.setArea;
setArea = strrep(setArea, 'AM', 'pAM');
setArea = strrep(setArea, 'ApAM', 'aAM');

%%
setK = 2:20; %paramClustering_global.setK; %Clustering.setK;
matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
    curK = setK(iK);
    matWSS(:,iK) = sum(Clustering_brainmask.resultKMeans(iK).SU_sumD); %sum(Clustering.resultKMeans(iK).SU_sumD);
    matWSS_roi(:,iK) = sum(Clustering_meanROI.resultKMeans(iK).roi_sumD); % for grouping ROIs for visualization purposes
end
totalSS = Clustering_brainmask.totalSS_SU;
propExplained = (totalSS-matWSS)./totalSS; %matExpVar./totalSS;
totalSS_roi = Clustering_meanROI.totalSS_roi;
propExplained_roi = (totalSS_roi-matWSS_roi)./totalSS_roi;
% Get the ROI sorting
curK_roi = 9;
locMode_roi = find(propExplained_roi(:,curK_roi-1)==mode(propExplained_roi(:,curK_roi-1)));
locMin_roi = find(propExplained_roi(:,curK_roi-1)==min(propExplained_roi(:,curK_roi-1)));
[sortedClust_roi, indSortROI] = sort(Clustering_meanROI.resultKMeans(curK_roi-1).roi_indCluster(:, locMode_roi(1)));

% Get the clustering results from voxel-based clustering
curK = 6; %8; %6; %10; %13; % 9; %6; %7;
locMode = find(propExplained(:,curK-1)==mode(propExplained(:,curK-1)));
locMin = find(propExplained(:,curK-1)==min(propExplained(:,curK-1)));
[sortedClust, indSortChan] = sort(Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(:, locMode(1)));
numROI = length(Clustering_meanROI.nameROI);
% orderROI = [1 2 22 3 4 35 34 12 13 14 29 30 6 7 8 36 9 10 11 32 15 5 23 26 37 27 28 16 17 33 18 19 20 21 31 24 25]; % 1:37;
% ordering of cells: ML-AF-AM-AAM, while maintaining the grouping
clear tempS
for iK = 1:curK
    tempS(iK).indSortChan_org = indSortChan(sortedClust==iK);
    tempS(iK).sortedClust_org = sortedClust(sortedClust==iK);
    tempS(iK).numCells = sum(sortedClust==iK);
    locML = find(Clustering_brainmask.infoCells.catAreaID(tempS(iK).indSortChan_org)==4);
    if isempty(locML)
        tempS(iK).indSortChan_reorder_MLfirst = tempS(iK).indSortChan_org;
    else
        tempS(iK).indSortChan_reorder_MLfirst = cat(1, tempS(iK).indSortChan_org(locML), tempS(iK).indSortChan_org(1:locML(1)-1));
    end
end
[~,reorderCluster] = sort(cat(1, tempS.numCells), 'descend');
% reorderCluster = 1:curK; %[8 1 9 5 3 7 4 10 6 2]; % ORder cell groups based on the number of neurons in each group (from largest to smallest)
indSortChan_reorder = cat(1, tempS(reorderCluster).indSortChan_org);
indSortChan_reorder_MLfirst = cat(1, tempS(reorderCluster).indSortChan_reorder_MLfirst);
sortedClust_reorder = cat(1, tempS(reorderCluster).sortedClust_org);
for iArea = 1:4
    ttt(iArea).indChanArea = find((Clustering_brainmask.infoCells.catAreaID==iArea)>0);
    ttt(iArea).indChanArea_rand = ttt(iArea).indChanArea(randperm(length(ttt(iArea).indChanArea)));
end
indChanArea_rand = cat(1,  ttt([2 1 3 4]).indChanArea_rand);
areaID_rand = Clustering_brainmask.infoCells.catAreaID(indChanArea_rand);

%% Curate the data for sharing
matR_neuron_fROI = matR(indChanArea_rand, indSortROI); % Fig 3A matrix

nameROI = Clustering_meanROI.nameROI(indSortROI);

infoCell.setArea = setArea(catAreaID(indChanArea_rand))'; % Fig 3A matrix: row information (face patch)

clusterID_org = Clustering_brainmask.resultKMeans(curK-1).SU_indCluster(indChanArea_rand, locMode(1));
[ta, tb] = sort(reorderCluster, 'ascend');
infoCell.groupID = tb(clusterID_org); % ID for 

ffsi = fsi.matFSI(locFaceCell, 1);

infoCell.fsi = ffsi(indChanArea_rand);

info.cell = infoCell;
info.ROI = nameROI;

save(fullfile(directory.dataHome, 'correlationMatrix_neuron_fROI.mat'), 'matR_neuron_fROI');
save(fullfile(directory.dataHome, 'info.mat'), 'info');

