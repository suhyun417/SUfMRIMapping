% computeCorrMap_AllCells_FPAreaSummary.m
%
% 2020/05/04 SHP
% Load the correlation matrix (unmasked, not pca-res) for all 368 cells
% from cortical face patches and compute 1) maximum of absolute correlation
% across cells in each area and 2) average map across cells in each area

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
nameSubjBOLD ='Art'; 
dirDataHome = fullfile(dirProcdata, 'parksh/_macaque');
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% load the masks (movie-driven & brain-only mask)
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp', 'brainMask_BlockAna3D');



% setNameSubjNeural{1} = {'Tor', 'Rho', 'Sig', 'Spi', 'Moc'}; % for AF: 11 12 13 14 15
% setNameSubjNeural{2} = {'Mat', 'Was'}; % for AM: 21 22
% setNameSubjNeural{3} = {'Dan', 'Moc'}; % for AM+: 31 32

%% Load the corr map and select valid channel: subject by subject
numSubject = size(setNameSubjNeural, 2);
matR_SU_all = [];
setArea = {'AF', 'AM', 'AM+', 'ML'};
for iSubj = 1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');
    
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
        case 'dan' %AM+
            subjID = 31;
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

corrMap_Area = struct([]);
for iArea = 1:length(unique(catAreaID))
    setR = matR_SU_all(:, catAreaID == iArea);
    matR_maxAbs = max(abs(setR), [], 2);
    matR_avg = mean(setR, 2);    
    
    corrMap_Area(iArea).nameArea = setArea{iArea};
    corrMap_Area(iArea).setSubjID = unique(catSubjID(catAreaID == iArea));
    corrMap_Area(iArea).setSubjName = unique(catSubjName(catAreaID == iArea, :));
    corrMap_Area(iArea).setChanID = catChanID(catAreaID == iArea);
    corrMap_Area(iArea).matR_max = matR_maxAbs;
    corrMap_Area(iArea).matR_avg = matR_avg;
end
    
corrMap_merged.catChanID = catChanID;
corrMap_merged.catSubjID = catSubjID;
corrMap_merged.catSubjID = catSubjName;
corrMap_merged.catAreaID = catAreaID;
corrMap_merged.matR = matR_SU_all;

%% save the file
save('/procdata/parksh/_macaque/CorrMap_SU_AllCellsArt_corticalFPMerged.mat', 'info*', 'corrMap_Area', 'corrMap_merged')
    
    
    
    
    


