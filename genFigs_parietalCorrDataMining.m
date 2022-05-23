% genFigs_parietalCorrDataMining.m
% 2022/04/18 SHP
% 1st playground to look at the face cells' correlation with voxels in the
% parietal cortex (a.k.a. "parietal project" with Marta Andujar)

% - read the excel file (/nifvault/projects/parksh/NeuroMRI/analysis/D99_Marta.xlsx)
% to retrieve which areas in D99 atlas we want to look at (Marta marked the
% areas)
% - get the coordinates from Artemis' D99 parcellations
% - load the face patch correlation map in Artemis brain

clear all;

%% Settings
flagBiowulf = 0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        directory.projects = '/Volumes/NIFVAULT/PROJECTS';
        directory.procdata = '/Volumes/NIFVAULT/PROCDATA';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/Volumes/NIFVAULT/LIBRARY';
        addpath(fullfile(directory.library, 'matlab_utils'));
    else % on virtual machine
        directory.projects = '/nifvault/projects';
        directory.procdata = '/nifvault/procdata';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/nifvault/library';
        addpath(fullfile(directory.library, 'matlab_utils'));
    end
end

%%
% get the ROIs to include from the spreadsheet
opts = spreadsheetImportOptions;
opts.VariableNames = {'indROI_obsolete', 'name', 'fullname', 'flagInclude', 'indROI'};
opts.VariableTypes = {'double', 'char', 'char', 'double', 'double'};
opts.DataRange = 'A2';

T = readtable(fullfile(directory.dataHome, 'D99_Marta.xlsx'), opts);
S = table2struct(T);
indParietalArea = S(cat(1, S.flagInclude)>0);

clear T

%% Get the voxel locations for each ROI
nameSubjBOLD = 'Art'; %'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataBOLD = fullfile(directory.dataHome, nameSubjBOLD);
load(fullfile(dirDataBOLD, sprintf('%sD99_ROIsegmentation.mat', nameSubjBOLD)), 'voltc_seg')

for iArea = 1:length(indParietalArea)
    locArea = find(voltc_seg==indParietalArea(iArea).indROI);
    indParietalArea(iArea).locVox = locArea;
end

% %% Any way we can do it in EPI space?
% for iArea = 1:length(indParietalArea)
%     clear locVox*
%     locVox_vol = voltc_seg==indParietalArea(iArea).indROI;
%     locVox_epi = permute(locVox_vol, [3 1 2]);
%     locVox_epi = decimate3D(locVox_epi, [3 3 3], .25);
%     
%     indParietalArea(iArea).locVox_epi = locVox_epi;
% end

%% checking the voxel locations that we are selecting (just for the sanity check)
setIndROI = cat(1, indParietalArea.indROI);
temp = ismember(voltc_seg, setIndROI); % voxels of any of the selected ROIs (UNION)

% Coronal
figure;
for iSlice = 85:140

    clf;
    ax1 = axes;
    im1 = imagesc(ax1, squeeze(voltc_seg(iSlice, :, :)));
    hold all;
    ax2 = axes;
    im2 = imagesc(ax2, squeeze(temp(iSlice, :, :)));
    im2.AlphaData = 0.7;
    linkaxes([ax1, ax2])
    ax2.Visible = 'off';
    input('')
end

% Sagittal
figure;
for iSlice = 2:100

    clf;
    ax1 = axes;
    im1 = imagesc(ax1, voltc_seg(:, :, iSlice));
    hold all;
    ax2 = axes;
    im2 = imagesc(ax2, temp(:, :, iSlice));
    im2.AlphaData = 0.7;
    linkaxes([ax1, ax2])
    ax2.Visible = 'off';
    input('')
end



%% Load correlation matrix of neurons 
setNameSubjNeural = {'Spi', 'Mat', 'Dan', 'Moc', 'Was', 'Dav'}; % Tor/Rho/Sig data cannot be converted to D99

%% Load the corr map and select valid channel: subject by subject
numSubject = size(setNameSubjNeural, 2);
% matR_SU_all = [];
% setArea = {'AF', 'AM', 'AAM', 'ML', 'NFP'};

% Prepare empty matrix for each area for concatenation across subjects
catCorrArea = cell(length(indParietalArea), 1); 

for iSubj = 1:numSubject
    
    clear matVolTc curMatCorr
    
    % load individual subject's correlation maps aligned to D99
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(directory.dataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new_D99resample.mat', nameSubjNeural, nameSubjBOLD)),...
        'catVoltc_new')
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)),...
        'paramCorr')
    
    % For Dango, include only face patch cells
    if strcmpi(nameSubjNeural, 'dan') % Dango has both AM+ and non-face patch cells, so assign different IDs for cells from each area
        locFP = cellfun(@isempty, strfind(cellstr(paramCorr.validChanID), 'NFP'));
        catVoltc_new = catVoltc_new(:,:,:,locFP);
    end
    
    matVolTc = reshape(catVoltc_new, prod(size(catVoltc_new, [1 2 3])), size(catVoltc_new, 4));
        
    for iArea = 1:length(indParietalArea)          
        curMatCorr = matVolTc(indParietalArea(iArea).locVox, :);
        catCorrArea{iArea} = cat(2, catCorrArea{iArea}, curMatCorr);
    end
    
    
    %% Subject & cell information related stuff: should be in a decent format to collate later
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
    % Include only face patch cells from Dango's
    if strcmpi(nameSubjNeural, 'dan') % Dango has both AM+ and non-face patch cells, so assign different IDs for cells from each area
        locFP = cellfun(@isempty, strfind(cellstr(paramCorr.validChanID), 'NFP'));
        infoPopulation_subj(iSubj).validChanIndex = [1:sum(locFP)]';
        infoPopulation_subj(iSubj).validChanID = paramCorr.validChanID(locFP, :);
        infoPopulation_subj(iSubj).validChan_subjID = ones(size(infoPopulation_subj(iSubj).validChanIndex)).*subjID(1);
    end
    %
    
end
    
% catCorr{iArea} = matCorr;
% catCorr_mean{iArea} = mean(matCorr);

% load('/nifvault/procdata/parksh/_macaque/tempParietalCorr.mat')


%% getting the cell info together across subjects
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

infoCells.catChanID = catChanID;
infoCells.catSubjName = catSubjName;
infoCells.catSubjID = catSubjID;
infoCells.catAreaID = catAreaID;
infoCells.setAreaID = setAreaID;

%%

matSCCorr_mean = cell2mat(cat(1, catCorr_mean(:)));
[catAreaName{1:length(indParietalArea)}] = deal(indParietalArea.name);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(matSCCorr_mean)
set(gca, 'CLim', [-1 1].*0.4)
colorbar
ylabel('ROI')
xlabel('Cells')
nCell_area = [12 40 48 30 12 21];
set(gca, 'XTick', cumsum(nCell_area))
set(gca, 'YTick', 1:length(indParietalArea), 'YTickLabel', catAreaName)





% load(fullfile(directory.dataHome, 'Art/Clustering_CorrMap_4FPs_faceselective_Movie123_probability.mat'), 'Clustering_brainmask', 'param*') 
% 
% % get only face selective cells
% load('/procdata/parksh/_macaque/multipleFP_fsi.mat')
% locFaceCell =  find(fsi.matFSI(:,1)>0.33); % find(abs(fsi.matFSI(:,1))>0.33);
% 
% matR = Clustering_meanROI.matR(locFaceCell, :);

% 
% indSubcorticalArea(1).name = 'PI';
% indSubcorticalArea(1).indROI = 196;
% indSubcorticalArea(2).name = 'PL';
% indSubcorticalArea(2).indROI = 197;
% indSubcorticalArea(3).name = 'PM';
% indSubcorticalArea(3).indROI = 198;
% indSubcorticalArea(4).name = 'Bi';
% indSubcorticalArea(4).indROI = 42;
% indSubcorticalArea(5).name = 'Bmc';
% indSubcorticalArea(5).indROI = 147;
% indSubcorticalArea(6).name = 'Bpc';
% indSubcorticalArea(6).indROI = 168;
% indSubcorticalArea(7).name = 'Ld';
% indSubcorticalArea(7).indROI = 9;
% indSubcorticalArea(8).name = 'Lv';
% indSubcorticalArea(8).indROI = 21;
% indSubcorticalArea(9).name = 'LGNm';
% indSubcorticalArea(9).indROI = 200;
% indSubcorticalArea(10).name = 'LGNp';
% indSubcorticalArea(10).indROI = 201;
% indSubcorticalArea(11).name = 'Claustrum';
% indSubcorticalArea(11).indROI = 176;
% indSubcorticalArea(12).name = 'Striatum';
% indSubcorticalArea(12).indROI = 175;

