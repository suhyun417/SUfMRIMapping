% testTemporalEvent.m
%
% Starting from a ROI, 1) sort cells based on the correlation values with
% that particular ROI (voxel), 2) see whether any specific temporal event
% determine the correlation level

S_neuralRegressor % do this first to get necessary structs such as STDPATH, dataBOLD, etc.

global dataBOLD S visROIs faceROIs DSP STDPATH 

% dirROI = fullfile(STDPATH.dataBOLD, 'ROIs');
% d_vis = dir(fullfile(dirROI, '*VisROIs.mat'));
% d_face = dir(fullfile(dirROI, '*faceROIs2.mat'));
% load(fullfile(dirROI, d_vis.name));
% load(fullfile(dirROI, d_face.name));
% 
% % get ROI information into structures with different names
% tempName_visROIs = char(fieldnames(load(fullfile(dirROI, d_vis.name))));
% tempName_faceROIs = char(fieldnames(load(fullfile(dirROI, d_face.name))));
% eval(['visROIs=', tempName_visROIs, ';']) % get visual ROI data into "visROIs"
% eval(['faceROIs=', tempName_faceROIs, ';']) % get face ROI data into "faceROIs"

ROIs = cat(2, visROIs', faceROIs'); % ROIs: 8x2 struct, 1st column visual ROIs, 2nd column face ROIs

% clear up
% eval(['clear ' tempName_visROIs])
% eval(['clear ' tempName_faceROIs])
% clear visROIs faceROIs


% Get ROI names & coordinates in EPI coords 
resizeFactor = DSP.proc.params3d.res./ROIs(1,1).params.res; % calculate the scaling factor

dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';


% Seed voxel correlation map 

% Get the averaged time series (across runs)
selMov=[1 2 3]; %1:9; %[1 2 3]; %4; %3; %2;

matBOLD = [];
matBOLD_ROI=[];

for iM = 1:length(selMov)
    idMov = selMov(iM);
    tLoc = find(dataBOLD.unimov==idMov);
    curvoltc = dataBOLD.mvoltc{tLoc};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    % 	subvoltc = curvoltc-repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    % 	pcvoltc = subvoltc./
    matBOLD = cat(4,matBOLD,pcvoltc);
end


flagFaceArea=1;
iFR=5; %4; %1:PO 2:PA 3:AM 4:AF 5:AL 6:MF 7:ML 8:PL

% for iFR=1:8 
clear voxROI a b c matR vecValidR

nameROI = ROIs(iFR,1+flagFaceArea).name(strfind(ROIs(iFR,1+flagFaceArea).name, '_')+1:end); %faceROIs(iFR).name(strfind(faceROIs(iFR).name, '_')+1:end);
voxROI=decimate3D(ROIs(iFR,1+flagFaceArea).vol3D, resizeFactor, .25); %decimate3D(faceROIs(iFR).vol3D, resizeFactor, .25); % turn anat_res ROIs into func_res ROIs
[a, b, c] = ind2sub(size(voxROI), find(voxROI==1)); % Get indices of AF voxels in EPI 3D coords
indVox_ROI = [a b c];
indVox_ROI_sub = sub2ind(size(voxROI), a, b, c);

% Find the highest correlated/involved voxel in one ROI
nameSubjNeural = 'Tor'; % 'Sig'; %'Tor';
nameSubjBOLD = 'Art'; %'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataNeural  = STDPATH.dataNeural;
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr')
% compute involvement/conjuction map
absMatR = abs(matR_SU);
absMatR_avg = mean(absMatR, 2);
% compute grand average of maps
matR_SU_grandAvg = mean(matR_SU, 2);

% quick check on the correlation distribution in this ROI
figure; hist(absMatR_avg(indVox_ROI_sub));
% get an index of voxel
[a ind]=sort(absMatR_avg(indVox_ROI_sub));
indVox = indVox_ROI_sub(ind(end));
[aaa bbb ccc] =  ind2sub(size(voxROI),indVox);
indVox_3d = [aaa bbb ccc];


% Distribution of correlation between all cells and this particular voxel
% (in this ROI)
vectorR_indVox = matR_SU(indVox, :);
[aa indCellSorted] = sort(vectorR_indVox, 'descend');

% crit = 0.1;
% indCellSorted_valid = indCellSorted(abs(aa)>crit);

indMovieNeuron = [1 2 3];

TR=2.4;
k = gampdf([-40:TR:40],4,2);

figTop5 = figure;
set(figTop5, 'Color', 'w', 'PaperPositionMode', 'auto')

for iChan = 1:3

    neuralrgrs=[];
    neuralrgrs = cat(1, S(paramCorr.validChanIndex(indCellSorted(iChan)), indMovieNeuron).mnFR);
    neuralrgrs = neuralrgrs-mean(neuralrgrs); % centering
    neuralrgrs = doConv(neuralrgrs,k);
    
    tempBOLD = squeeze(matBOLD(indVox_3d(1), indVox_3d(2), indVox_3d(3), :));
    
    tsBOLD_norm = (tempBOLD-nanmean(tempBOLD))./nanstd(tempBOLD);
    tsNeural_norm = zscore(neuralrgrs);
    
    figure(figTop5)
    subplot(3, 1, iChan);
    plot(tsBOLD_norm.*(-1), 'k-', 'LineWidth', 2)
    hold on
    plot(tsNeural_norm, 'm-', 'LineWidth', 2)
    title(sprintf('Cell ID: %s, \rho = %2.2f', paramCorr.validChanID(indCellSorted(iChan), :), aa(iChan)))
    
end

figBottom5 = figure;
set(figBottom5, 'Color', 'w', 'PaperPositionMode', 'auto')

indCellSorted_flipped = fliplr(indCellSorted);
corr_flipped = fliplr(aa);
for iChan = 1:3
    
    neuralrgrs=[];
    neuralrgrs = cat(1, S(paramCorr.validChanIndex(indCellSorted_flipped(iChan)), indMovieNeuron).mnFR);
    neuralrgrs = neuralrgrs-mean(neuralrgrs); % centering
    neuralrgrs = doConv(neuralrgrs,k);
    
    tempBOLD = squeeze(matBOLD(indVox_3d(1), indVox_3d(2), indVox_3d(3), :));
    
    tsBOLD_norm = (tempBOLD-nanmean(tempBOLD))./nanstd(tempBOLD);
    tsNeural_norm = zscore(neuralrgrs);
    
    figure(figBottom5)
    subplot(3, 1, iChan);
    plot(tsBOLD_norm.*(-1), 'k-', 'LineWidth', 2)
    hold on
    plot(tsNeural_norm, 'm-', 'LineWidth', 2)
    title(sprintf('Cell ID: %s, \rho = %2.2f', paramCorr.validChanID(indCellSorted_flipped(iChan), :), corr_flipped(iChan)))
    
end




