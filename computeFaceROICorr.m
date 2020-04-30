
clear all;


nameSubjNeural = 'Tor'; % 'Sig'; %'Tor';
nameSubjBOLD = 'Art'; %'Ava'; %'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';

%
addpath('/library/matlab_utils/')

% Load files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% 1) fMRI correlation maps
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr')

% 2) Clustering results
% load(fullfile(dirDataNeural, 'Clustering_TorArtMovie123_new.mat')) %
load(fullfile(dirDataNeural, 'Clustering_TorArtMovie123_new_masked.mat')) % clustering based on the maps with movie-driven mask applied
% load(fullfile(dirDataNeural, 'Clustering_TorArtAvaMovie123.mat')) % based on two sets of whole brain maps from two different monkeys
% Clustering = ClusteringAll;

% 3) fMRI movie-driven activity mask
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');

% 4) Single-unit time courses
load(fullfile(dirDataNeural, sprintf('%s_movieTS_SU_indMov.mat', nameSubjNeural)))

% filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% fileNameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
% filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];
%
% fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
% load(fullfile(dirDataNeural, filenameNeural))
% load(fullfile(dirDataNeural, fileNameNeural_BLP))
% fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
% load(fullfile(dirDataBOLD, filenameBOLD))

% % Get movie IDs common in two dataset
% commonSetMovie = intersect(paramBOLD.unimov, paramSDF.setMovIDs);
% dataBOLD.mvoltc = voltcIndMov(commonSetMovie);
% dataBOLD.unimov = commonSetMovie;
%
% setMovie = [1 2 3];

%% 1. Average fMRI maps for each cluster
matIndClust_SU = cat(2, Clustering_moviemask.resultKMeans.SU_indCluster);
% matIndClust_SU = cat(2, Clustering.resultKMeans.SU_indCluster);

sortTargetK = 5; %6:8 %4; %:5 %sortTargetK = 4; %8; %6; %7; %6; %7;
[sortedClust, indSortChan]=sort(matIndClust_SU(:,sortTargetK-1));

% average maps for each cluster
nx = 40; ny = 64; nz = 32;
mapR_Cluster=[]; R_avgCluster=[];

for iK = 1:sortTargetK
    R_avgCluster(:,iK) = mean(matR_SU(:,indSortChan(sortedClust==iK)), 2);
    tempClustMapR = reshape(R_avgCluster(:,iK), [nx, ny, nz]);
    mapR_Cluster = cat(4, mapR_Cluster, tempClustMapR);
end


%% 2. ROIs
dirROI = fullfile(dirDataBOLD, 'ROIs');
% if ~exist(dirROI)
%     visROIs=[]; faceROIs=[];
% else
load(fullfile(dirROI, 'e66_faceROIs_L.mat'))
load(fullfile(dirROI, 'e66_faceROIs2.mat'))

faceROIs = struct([]);
faceROIs = cat(2, e66_faceROIs2, e66_faceROIs_L);

clear e66*

% Get ROI names & coordinates in EPI coords 
resizeFactor = [3 3 3]; % DSP.proc.params3d.res./ROIs(1,1).params.res; % calculate the scaling factor

% switch lower(nameSubjBOLD)
%     case 'art'
%         indAF = 4; 
%     case 'ava'
%         indAF = 7;
% end
setCluster = [3 4 5 2 1];
meanCorr_ROI=[]; steCorr_ROI=[];
for iR = 1:length(faceROIs) 
nameROI = faceROIs(iR).name; %(strfind(faceROIs(iR).name, '_')+1:end); %faceROIs(iFR).name(strfind(faceROIs(iFR).name, '_')+1:end);
voxROI=decimate3D(faceROIs(iR).vol3D, resizeFactor, .25); %decimate3D(faceROIs(iFR).vol3D, resizeFactor, .25); % turn anat_res ROIs into func_res ROIs
[a, b, c] = ind2sub(size(voxROI), find(voxROI==1)); % Get indices of AF voxels in EPI 3D coords
indVox_ROI = [a b c];
indVox_ROI_sub = sub2ind(size(voxROI), a, b, c);

matCorr_ROI=[];
for iVox = 1:length(a)
    matCorr_ROI(:,iVox) = mapR_Cluster(a(iVox), b(iVox), c(iVox), :); %matBOLD_shuffle(a(iVox), b(iVox), c(iVox),:); %matBOLD(a(iVox), b(iVox), c(iVox),:);
end
meanCorr_ROI(:,iR) = nanmean(matCorr_ROI, 2);
steCorr_ROI(:,iR) = nanstd(matCorr_ROI, [], 2)./sqrt(size(matCorr_ROI,2)-1);

figure(100)
subplot(length(faceROIs), 1, iR)
plot(meanCorr_ROI(setCluster, iR), 'bo')
hold on;
line([1:5; 1:5], [meanCorr_ROI(setCluster, iR) - steCorr_ROI(setCluster,iR), meanCorr_ROI(setCluster, iR) + steCorr_ROI(setCluster,iR)]', 'Color', 'b')
end


setROIs_right = [3 4 5 6 7 8];
setROIs_left = [10 11 12 14 13 15];

ClusterCorr_faceROI.nameROI = cat(1, faceROIs.name);
ClusterCorr_faceROI.meanCorr_ROI = meanCorr_ROI;
ClusterCorr_faceROI.steCorr_ROI = steCorr_ROI;
ClusterCorr_faceROI.setROIs_right = setROIs_right;
ClusterCorr_faceROI.setROIs_left = setROIs_left;
ClusterCorr_faceROI.indNewClusterOrder = setCluster;

save(fullfile(dirDataNeural, 'ClusterCorrMap_FaceROI.mat'), 'ClusterCorr_faceROI')


