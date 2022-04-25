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
opts.VariableNames = {'indROI', 'name', 'fullname', 'flagInclude'};
opts.VariableTypes = {'double', 'char', 'char', 'double'};
opts.DataRange = 'A2';

T = readtable(fullfile(directory.projects, '/parksh/NeuroMRI/analysis/D99_Marta.xlsx'), opts);
S = table2struct(T);
indParietalArea = S(cat(1, S.flagInclude)>0);

clear T

%% Get all the 
nameSubjBOLD = 'Art'; %'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataBOLD = fullfile(directory.dataHome, nameSubjBOLD);
load(fullfile(dirDataBOLD, sprintf('%sD99_ROIsegmentation.mat', nameSubjBOLD)), 'voltc_seg')
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


for iArea = 1:length(indParietalArea)
    locArea = find(voltc_seg==indParietalArea(iArea).indROI);
    indParietalArea(iArea).locVox = locArea;
end

for iArea = 1:11 %length(indSubcorticalArea)
    matCorr = [];
    for iSubj = 2:6
        nameSubjNeural = setNameSubjNeural{iSubj};
        dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
        
        load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie%s_new_D99resample.mat', nameSubjNeural, nameSubjBOLD, MovieStr)),...
            'catVoltc_new')
        
        nCell = size(catVoltc_new,4);
        
        
        for iCell=1:nCell
            curVoltc = catVoltc_new(:,:,:,iCell);
            matCorr = cat(2, matCorr, curVoltc(indParietalArea(iArea).locVox));
        end
        
    end
    
    % add pulvinar neurons
    nameSubjNeural = 'Dex';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
        
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie1_new_D99resample.mat', nameSubjNeural, nameSubjBOLD)),...
        'catVoltc_new')
    
    nCell = size(catVoltc_new,4);
        
    for iCell=1:nCell
        curVoltc = catVoltc_new(:,:,:,iCell);
        matCorr = cat(2, matCorr, curVoltc(indParietalArea(iArea).locVox));
    end
    
    %
    catCorr{iArea} = matCorr;
    catCorr_mean{iArea} = mean(matCorr);
end

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
