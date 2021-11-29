% genFig_multipleFP_Revision_faceRegressorCellVoxels.m
%
% 2021/11/16 SHP
% Generate a figure showing face regressor (time course of face appearance
% during the movie) with example neighboring neurons and voxels response
% time course, to demonstrate responses and correlations that we are
% dealing with 
%   - This is a response to a Reviewer's point during revision from Science
%   Advances

clear all;

%% Settings
flagBiowulf = 1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
%     addpath('/data/parks20/analysis/NeuroMRI/'); % to use doConv.m function
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        directory.projects = '/Volumes/NIFVAULT/projects';
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
% dirFig = fullfile(directory.projects, 'parksh/NeuroMRI/_labNote/_figs');

load(fullfile(directory.dataHome, 'tempRevisionAnalysis_faceRegressorCellVoxel.mat'))

flagBiowulf = 1;
if flagBiowulf
clear directory
directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
end


% %% 1. Retrieve all the necessary data
% %% 1-1. Load face regressor
% % dataDir ='/procdata/parksh/MovieRegressors/'; %'/Volumes/PROCDATA/parksh/MovieRegressors/'; %'/procdata/parksh/MovieRegressors/';
% % load(fullfile(directory.procdata, 'parksh/MovieRegressors/HighLevelRGR.mat'))
% 
% flagSM = 1; % flag for compression and smoothing
% setMovie = [1 2 3];
% fullRGR4fps = createMovieRGR_4fps_indMov(setMovie, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
% 
% tempCat = cat(1, fullRGR4fps(:).regressors);
% % faceRegressor = tempCat(:, 20); % all the moments of presentation of faces in full and side views
% faceRegressor = resample(tempCat(:, 20), 0.25*100, 2.4*100); %matRGR = resample(catRGR, 0.25*100, 2.4*100);
% 
% % MION function
% TR=2.4;
% k = gampdf([-40:TR:40],4,2);
% faceRegressor_mion = doConv(faceRegressor, k); % 
% 
% 
%% 1-2.Load cell time series
load(fullfile(directory.dataHome, 'matSDF_Movie123_allCells.mat'), 'matTS_FP')

setCellIDs = {'06Dav', '25Dav', '27Dav'};
tLoc = find(contains(matTS_FP.catChanID, setCellIDs));  

cellTS = matTS_FP.matNeuralRGR(:, tLoc);
cellTS_norm = (cellTS-nanmean(cellTS))./nanstd(cellTS);
% 
% %% 1-3.Load voxel time series
% % 1-3-1. Load entire brain
% load(fullfile(directory.dataHome, 'Art', 'Art_movieTS_fMRI_indMov.mat')); % '_movieTS_fMRI_indMov.mat'];
% 
% dataBOLD.mvoltc = voltcIndMov([1 2 3]);
% 
% % 1. fMRI tc
% fmritc=[];
% for iM = 1:3
%     curvoltc = dataBOLD.mvoltc{iM};
%     avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
%     if ~isempty(find(avgvoltc==0, 1))
%         avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
%     end
%     pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
%     fmritc = cat(4,fmritc,pcvoltc);
% end

% 2. ROIs
% nameSubjBOLD = 'Art';
% 
% % load fmri ROIs
% dirDataHome = '/procdata/parksh/';
% dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
% dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

dirROI = fullfile(directory.dataHome, 'Art/ROIs');

d_vis = dir(fullfile(dirROI, '*VisROIs.mat'));
d_face = dir(fullfile(dirROI, '*faceROIs2.mat'));
load(fullfile(dirROI, d_vis.name));
load(fullfile(dirROI, d_face.name));

% get ROI information into structures with different names
tempName_visROIs = char(fieldnames(load(fullfile(dirROI, d_vis.name))));
tempName_faceROIs = char(fieldnames(load(fullfile(dirROI, d_face.name))));
eval(['visROIs=', tempName_visROIs, ';']) % get visual ROI data into "visROIs"
eval(['faceROIs=', tempName_faceROIs, ';']) % get face ROI data into "faceROIs"

% clear up
eval(['clear ' tempName_visROIs])
eval(['clear ' tempName_faceROIs])

% Get ROI names & coordinates in EPI coords 
setIndROI = [1 1; 1 3; 2 7]; % visROIs(1): V1fovea, visROIs(3): MT, faceROIs(7): face patch ML

[dvolx, dvoly, dvolz]  = size(visROIs(1).vol3D); % dimension of volume
resizeFactor = [3 3 3]; % DSP.proc.params3d.res./faceROIs(1).params.res; % calculate the scaling factor

for iROI = 1:size(setIndROI, 1)
    curROI = struct([]);
    switch setIndROI(iROI, 1)
        case 1
            curROI = visROIs(setIndROI(iROI, 2));
        case 2
            curROI = faceROIs(setIndROI(iROI, 2));
    end
    
    nameROI = curROI.name;
    voxROI = decimate3D(curROI.vol3D, resizeFactor, .25);
    [a, b, c] = ind2sub(size(voxROI), find(voxROI==1)); % Get indices of AF voxels in EPI 3D coords
    indVox_ROI = [a b c];
    indVox_ROI_sub = sub2ind(size(voxROI), a, b, c);
        
    matTS_ROI=[];
    for iVox = 1:length(a)
        matTS_ROI(:,iVox) = fmritc(a(iVox), b(iVox), c(iVox), :); %matBOLD_shuffle(a(iVox), b(iVox), c(iVox),:); %matBOLD(a(iVox), b(iVox), c(iVox),:);
    end
    
    matTS(iROI).matTS_ROI = matTS_ROI.*(-1);
    matTS(iROI).meanTS_ROI = mean(matTS_ROI, 2);
end



%% Plot
% faceRegressor, faceRegressor_mion
% cellTS, cellTS_norm
% matTS(iROI).matTS_ROI(:, 10), matTS(iROI).meanTS_ROI
% fmriTS = cat(2, matTS.meanTS_ROI); %cat(2, matTS(1).matTS_ROI(:,10), matTS(2).matTS_ROI(:,10), matTS(3).matTS_ROI(:, 7));
fmriTS = cat(2, matTS(1).matTS_ROI(:,7), matTS(2).matTS_ROI(:,45), matTS(3).matTS_ROI(:, 20));

fmriTS_norm = (fmriTS-nanmean(fmriTS))./nanstd(fmriTS);

cMap_Area = [179 226 205; 141 160 203; 252 141 98; 231 41 138]./255; % AF-pAM-aAM-ML

%% example 3 cells with similar time course 
taxis = 2.4:2.4:900;
cMap_fmri = cool(3); %ones(3).*([0.77 0.5 0]'); %winter(3); %cool(3); %ones(3).*([0.7 0.3 0]'); % gray scale for three fMRI voxels from different ROIs
for iCell = 1:3
    figure;
    set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 100 650 160]);
    pp = plot(taxis, fmriTS_norm, 'LineWidth', 1);%
    hold on;
    plot(taxis, cellTS_norm(:, iCell), 'k-', 'LineWidth', 1);
    set(pp, {'Color'}, mat2cell(cMap_fmri, [1 1 1], [3]))   
    ylim([-4 4.5])
    
    set(gca, 'TickDir', 'out', 'box', 'off', 'XColor', 'k', 'YColor', 'k')
    print(gcf, fullfile(dirFig, sprintf('multipleFP_sFig_cellVoxTS_%s', setCellIDs{iCell})), '-depsc')
end

%% Correlation values
matR = corr(cellTS, fmriTS, 'rows', 'complete', 'type', 'spearman');
% matR = matR*(-1); %becuase of the MION

%% face regressor ts
taxis = 2.4:2.4:900;

faceRegressor_mion([1:7,126:126+7, 251:251+7]) = NaN;
faceRegressor_mion_norm = (faceRegressor_mion-nanmean(faceRegressor_mion))./nanstd(faceRegressor_mion);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [1200 100 650 160]);
plot(taxis, faceRegressor_mion_norm, '-', 'LineWidth', 1, 'Color', [154 205 50]./255);
set(gca, 'TickDir', 'out', 'box', 'off', 'XColor', 'k', 'YColor', 'k')
set(gca, 'YColor', 'none')
print(gcf, fullfile(dirFig, 'multipleFP_sFig_faceRegressor_mion'), '-depsc')

