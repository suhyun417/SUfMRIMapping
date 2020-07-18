% computeCorrMap_AllCells_FPAreaSummary_pcares.m
%
% 2020/07/09 SHP: compute grand summary of all face patch areas (mean and
% sum of individual summary maps)
% 2020/06/29 SHP: include new ML neurons and Dango's peri-AM neurons
% 2020/05/04 SHP
% Load the correlation matrix (unmasked, not pca-res) for all 368 cells
% from cortical face patches and compute 1) maximum of absolute correlation
% across cells in each area and 2) average map across cells in each area
% and 3) for each voxel, fraction of neurons that have the correlation
% higher than certain criterion (e.g. 0.3)

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
    
% Add necessary toolbox 
addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution

% Set directories 
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was', 'Dav'};
nameSubjBOLD = 'Art'; %'Ava'; %'Art';  %'Ava'; % 'Art'; 
dirDataHome = fullfile(dirProcdata, 'parksh/_macaque');
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% % load the masks (movie-driven & brain-only mask)
% load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');



% setNameSubjNeural{1} = {'Tor', 'Rho', 'Sig', 'Spi', 'Moc'}; % for AF: 11 12 13 14 15
% setNameSubjNeural{2} = {'Mat', 'Was'}; % for AM: 21 22
% setNameSubjNeural{3} = {'Dan', 'Moc'}; % for AM+: 31 32

%% Load the corr map and select valid channel: subject by subject
numSubject = size(setNameSubjNeural, 2);
matR_SU_all = [];
setArea = {'AF', 'AM', 'AAM', 'ML', 'NFP'};
for iSubj = 1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new_pcares.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');
    
    % For each cell, save the subject ID: Mochi's AF cells get 15, Mochi's
    % AM cells get 32
    switch lower(nameSubjNeural) % for AF: 11 12 13 14 15; for AM: 21 22; for AM+: 31 32; for ML: 41
        case 'tor'
            subjID = 11;
        case 'rho'
            subjID = 12;
        case 'sig'
            subjID = 13;
        case 'spi'
            subjID = 14;
        case 'mat' %AM
            subjID = 21;
        case 'dan' %AM+ and peri-AM+ (no face patch)
            subjID = [31 100];
        case 'moc' % AF and AM+
            subjID = [15 32];
        case 'was' %AM
            subjID = 22;
        case 'dav' %ML
            subjID = 41;
    end
        
    infoPopulation_subj(iSubj).nameSubj = nameSubjNeural;
    infoPopulation_subj(iSubj).validChanIndex = (1:length(paramCorr.validChanIndex))'; %validChanIndex_clustering;
    infoPopulation_subj(iSubj).validChanID = paramCorr.validChanID; %validChanID_clustering;
    infoPopulation_subj(iSubj).validChan_subjID = ones(size(infoPopulation_subj(iSubj).validChanIndex)).*subjID(1);
     if strcmpi(nameSubjNeural, 'moc') % Mochi has both AF and AM+ cells, so assign different IDs for cells from each area
        locAM = ~cellfun(@isempty, strfind(cellstr(paramCorr.validChanID), 'AM'));
        infoPopulation_subj(iSubj).validChan_subjID(locAM) = subjID(2);
     end
     if strcmpi(nameSubjNeural, 'dan') % Dango has both AM+ and non-face patch cells, so assign different IDs for cells from each area
        locNFP = ~cellfun(@isempty, strfind(cellstr(paramCorr.validChanID), 'NFP'));
        infoPopulation_subj(iSubj).validChan_subjID(locNFP) = subjID(2);
     end
    
    matR_SU_valid = matR_SU(:, infoPopulation_subj(iSubj).validChanIndex);
    matR_SU_all = cat(2, matR_SU_all, matR_SU_valid);
    
     clear matR_SU matR_SU_valid            
end

curTS = struct([]);
for iC = 1:length(infoPopulation_subj)
    curTS(iC).validChanID = strcat(cellstr(infoPopulation_subj(iC).validChanID), infoPopulation_subj(iC).nameSubj); %cellstr(paramClustering(iC).validChanID);
    curTS(iC).setSubjName = cellstr(repmat(infoPopulation_subj(iC).nameSubj, length(infoPopulation_subj(iC).validChanID), 1));
    %     indMonkey = cat(1, indMonkey, ones(size(paramClustering(iC).validChanIndex)).*iC);
end
catChanID = cat(1, curTS.validChanID);
catSubjName = cat(1, curTS.setSubjName);
catSubjID = cat(1, infoPopulation_subj.validChan_subjID);
catAreaID = floor(catSubjID./10);
setAreaID = unique(catAreaID);

corrMap_Area = struct([]);
for iArea = 1:length(setAreaID)
    idArea = setAreaID(iArea);
    setR = matR_SU_all(:, catAreaID == idArea);
    matR_maxAbs = max(abs(setR), [], 2);
    matR_avg = mean(setR, 2);    

    critCorr = 0.3;
    matValidVox = abs(setR)>critCorr;
    fractionHighCorrCell = sum(matValidVox, 2)./size(matValidVox,2);
    
    corrMap_Area(iArea).nameArea = setArea{iArea};
    corrMap_Area(iArea).setSubjID = unique(catSubjID(catAreaID == idArea));
    corrMap_Area(iArea).setSubjName = unique(catSubjName(catAreaID == idArea, :));
    corrMap_Area(iArea).setChanID = catChanID(catAreaID == idArea);
    corrMap_Area(iArea).catSubjID = catSubjID(catAreaID==idArea);
    corrMap_Area(iArea).matR = setR;
    corrMap_Area(iArea).matR_max = matR_maxAbs;
    corrMap_Area(iArea).matR_avg = matR_avg;
    corrMap_Area(iArea).critCorr = critCorr;
    corrMap_Area(iArea).fractionHighCorrCell = fractionHighCorrCell;
end


%% Merge all the neurons and all the areas together
corrMap_merged.catChanID = catChanID;
corrMap_merged.catSubjID = catSubjID;
corrMap_merged.catSubjName = catSubjName;
corrMap_merged.catAreaID = catAreaID;
corrMap_merged.matR = matR_SU_all;

% average fraction of neurons across areas
catMap_fraction = cat(2, corrMap_Area.fractionHighCorrCell);
corrMap_merged.meanFractionAcrossArea = mean(catMap_fraction, 2);

%% Merge only face patches
setAreaFP = 1:4;

corrMap_merged_FP.catChanID = cat(1, corrMap_Area(setAreaFP).setChanID); % catChanID;
corrMap_merged_FP.catSubjID = cat(1, corrMap_Area(setAreaFP).catSubjID); % catSubjID;
corrMap_merged_FP.setArea = {corrMap_Area(setAreaFP).nameArea}; % catSubjName;
corrMap_merged_FP.catAreaID = floor(corrMap_merged_FP.catSubjID./10); % catAreaID;
corrMap_merged_FP.matR = cat(2, corrMap_Area(setAreaFP).matR); % matR_SU_all;

% average fraction of neurons across areas
catMap_fraction_FP = cat(2, corrMap_Area(setAreaFP).fractionHighCorrCell);
catMap_matR_max_FP = cat(2, corrMap_Area(setAreaFP).matR_max);
catMap_matR_avg_FP = cat(2, corrMap_Area(setAreaFP).matR_avg);
corrMap_merged_FP.meanFractionAcrossArea = mean(catMap_fraction_FP, 2);
corrMap_merged_FP.meanMaxAbs = mean(catMap_matR_max_FP, 2);
corrMap_merged_FP.grandMeanR = mean(catMap_matR_avg_FP, 2);


%% save the file
save(sprintf('/procdata/parksh/_macaque/CorrMap_SU_AllCells%s_corticalFPMerged_pcares.mat', nameSubjBOLD), ...
    'info*', 'corrMap_*')
    
    
    
    
    


