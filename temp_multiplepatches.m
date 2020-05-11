% scribble things regarding multiple patch analyses until an idea gets
% shaped

clear all;

%%
addpath('/library/matlab_utils/')

%% Subcortical maps
ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS';
    dirProcdata = '/Volumes/PROCDATA';
    dirDataHome = fullfile(dirProcdata, 'parksh', '_macaque');
    dirLibrary = '/Volumes/LIBRARY';
    addpath(fullfile(dirLibrary, 'matlab_utils'));
else % on virtual machine
    dirProjects = '/projects';
    dirProcdata = '/procdata';
    dirDataHome = fullfile(dirProcdata, 'parksh', '_macaque');
    dirLibrary = '/library';
    addpath(fullfile(dirLibrary, 'matlab_utils'));
end

nameSubjBOLD = 'Art';
load(fullfile(dirDataHome, sprintf('%sD99_ROIsegmentation.mat', nameSubjBOLD)), 'voltc_seg')

indSubcorticalArea(1).name = 'PI';
indSubcorticalArea(1).indROI = 196;
indSubcorticalArea(2).name = 'PL';
indSubcorticalArea(2).indROI = 197;
indSubcorticalArea(3).name = 'PM';
indSubcorticalArea(3).indROI = 198;
indSubcorticalArea(4).name = 'Bi';
indSubcorticalArea(4).indROI = 42;
indSubcorticalArea(5).name = 'Bmc';
indSubcorticalArea(5).indROI = 147;
indSubcorticalArea(6).name = 'Bpc';
indSubcorticalArea(6).indROI = 168;
indSubcorticalArea(7).name = 'Ld';
indSubcorticalArea(7).indROI = 9;
indSubcorticalArea(8).name = 'Lv';
indSubcorticalArea(8).indROI = 21;
indSubcorticalArea(9).name = 'LGNm';
indSubcorticalArea(9).indROI = 200;
indSubcorticalArea(10).name = 'LGNp';
indSubcorticalArea(10).indROI = 201;
indSubcorticalArea(11).name = 'Claustrum';
indSubcorticalArea(11).indROI = 176;
indSubcorticalArea(12).name = 'Striatum';
indSubcorticalArea(12).indROI = 175;

for iArea = 1:length(indSubcorticalArea)
    locArea = find(voltc_seg==indSubcorticalArea(iArea).indROI);
    indSubcorticalArea(iArea).locVox = locArea;
end

load(fullfile(dirDataHome, 'CorrMap_SU_AllCellsArt_corticalFPMerged_D99resample.mat'),...
     'catVoltc_new', 'infoFile_org');
setFname = {infoFile_org.name}';
setFPArea = {'ML', 'AF', 'AM', 'AAM', 'MeanAcAreas'};

indFract = contains(setFname, 'FractHighCorrCell');
setFname_fract = setFname(indFract);

for iArea = 1:length(indSubcorticalArea)
    
    matFract = [];
    for iFPArea = 1:length(setFPArea)        
        idFPArea = find(contains(setFname_fract, setFPArea{iFPArea}));
        
        curVoltc = catVoltc_new(:,:,:,idFPArea);
        matFract = cat(2, matFract, curVoltc(indSubcorticalArea(iArea).locVox));
    end
        
    %
    catFract{iArea} = matFract;
    catFract_mean{iArea} = mean(matFract);
end

matSCFract = cell2mat(cat(1, catFract_mean(:)));
[catAreaName{1:length(indSubcorticalArea)}] = deal(indSubcorticalArea.name);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(matSCFract')
set(gca, 'CLim', [0 1].*0.1)
colorbar
xlabel('ROI')
ylabel('Areas')
set(gca, 'XTick', 1:length(indSubcorticalArea), 'XTickLabel', catAreaName)
set(gca, 'YTick', 1:length(setFPArea), 'YTickLabel', setFPArea) %cumsum(nCell_area))
title('Fraction of high correlation (r>0.3) neurons for each subcortical ROIs')
% print(gcf, fullfile(dirFig, 'PulvinarAmygdala_DavTorSpiMatDan_Movie123'), '-r150', '-dtiff')



