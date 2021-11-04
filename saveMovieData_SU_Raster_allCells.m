% saveMovieData_SU_Raster_allCells.m
%
% 2021/11/1 SHP
% Modified from saveMovieData_SU_SDF_allCells.m
% Prepare time-series matrices for each face patch across different animals for all the neurons
%       - rasters (1ms bin)


clear all;

%% Settings
flagBiowulf = 0;

if flagBiowulf
    directory.dirDataHome = '/data/parks20/procdata/NeuroMRI/';
    addpath('/data/parks20/analysis/NeuroMRI/'); % to use doConv.m function
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

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was', 'Dav'};

%% Parameters for time series
setMovie = [1 2 3]; 

setArea = {'AF', 'AM', 'AAM', 'ML', 'NFP'};

paramTS.setMovie = setMovie;

%% Load and compute rasters and SDFs for all cells
for iSubj = 1:length(setNameSubjNeural) %curSetNameSubjNeural)
    
    % Load the data
    nameSubjNeural = setNameSubjNeural{iSubj}; %curSetNameSubjNeural{iSubj};
    dirDataNeural = fullfile(directory.dataHome, nameSubjNeural);
    
    filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
    load(fullfile(dirDataNeural, filenameNeural))
    fprintf(1, '\nLoading single unit data of %s: %s ....\n ', nameSubjNeural, filenameNeural)
    
    % For each cell, save the subject ID that indicates area and subject: Mochi's AF cells get 15, Mochi's AM cells get 32
    switch lower(nameSubjNeural) % for AF: 11 12 13 14 15; for AM: 21 22; for AM+: 31 32
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
    
    [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
    validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
    
    if isempty(validC)
        continue;
    end  
    
    infoTS_subj(iSubj).nameSubj = nameSubjNeural;
    infoTS_subj(iSubj).validChanID = paramSDF.setCellIDs(validC);
    infoTS_subj(iSubj).validChanIndex_org = validC;
%     infoTS(iSubj).validChanIndex_cat = (1:length(validC))';
    infoTS_subj(iSubj).validChan_subjID = ones(size(infoTS_subj(iSubj).validChanID)).*subjID(1);
    if strcmpi(nameSubjNeural, 'moc') % Mochi has both AF and AM cells, so need to select cells from each area
        locAM = ~cellfun(@isempty, strfind(cellstr(infoTS_subj(iSubj).validChanID), 'AM'));
        infoTS_subj(iSubj).validChan_subjID(locAM) = subjID(2);
    end
    if strcmpi(nameSubjNeural, 'dan') % Dango has both AM+ and non-face patch cells, so assign different IDs for cells from each area
        locNFP = ~cellfun(@isempty, strfind(cellstr(infoTS_subj(iSubj).validChanID), 'NFP'));
        infoTS_subj(iSubj).validChan_subjID(locNFP) = subjID(2);
    end
    
%     if iSubj < 4 % previous recordings with no fingerprinting data
%         continue;
%     end
       
    indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
    % 1. Raster of spikes (spike count of 1ms bin)
    sizeTimeBin_sec = 0.001;
    raster = createCellRegressor_indMov_discreteTime(dirDataNeural, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
        setMovie, sizeTimeBin_sec); 
    % concatenate across movies
    matRaster=[];
    for iUnit = 1:size(raster,1)
        tempFR = cat(1, raster(iUnit, :).mnFR);
        matRaster(:,iUnit) = tempFR;
    end
    matTS_subj(iSubj).matRaster = matRaster;

    
end



%% Concatenate across subjects
ttt = struct([]);
for iC = 1:length(infoTS_subj)
    ttt(iC).validChanID = strcat(cellstr(infoTS_subj(iC).validChanID), infoTS_subj(iC).nameSubj); %cellstr(paramClustering(iC).validChanID);
    ttt(iC).setSubjName = cellstr(repmat(infoTS_subj(iC).nameSubj, length(infoTS_subj(iC).validChanID), 1));
    %     indMonkey = cat(1, indMonkey, ones(size(paramClustering(iC).validChanIndex)).*iC);
end
catChanID = cat(1, ttt.validChanID);
catSubjName = cat(1, ttt.setSubjName);
catSubjID = cat(1, infoTS_subj.validChan_subjID);
catAreaID = floor(catSubjID./10);
setAreaID = unique(catAreaID);

matTS_all.matRaster = cat(2, matTS_subj.matRaster);
% matTS_all.matSDF = cat(2, matTS_subj.matSDF);

matTS_all.infoCell.catChanID = catChanID;
matTS_all.infoCell.catSubjName = catSubjName;
matTS_all.infoCell.catSubjID = catSubjID;
matTS_all.infoCell.catAreaID = catAreaID;

% %% For each area (recording sites), concatenate across subject
% cat_matRaster = matTS_all.matRaster;
% cat_matNeuralRGR = matTS_all.matNeuralRGR;
% 
% % matTS_area = struct([]);
% setArea = {'AF', 'AM', 'AAM', 'ML', 'NFP'};
% for iArea = 1:length(setAreaID)
%     idArea = setAreaID(iArea);
%     
%     matTS_area(iArea).nameArea = setArea{iArea};
%     matTS_area(iArea).setSubjID = unique(catSubjID(catAreaID == idArea));
%     matTS_area(iArea).setSubjName = unique(catSubjName(catAreaID == idArea, :));
%     matTS_area(iArea).setChanID = catChanID(catAreaID == idArea);
%     matTS_area(iArea).catSubjID = catSubjID(catAreaID==idArea);
%     matTS_area(iArea).matFR_TR = cat_matRaster(:, catAreaID == idArea);
%     matTS_area(iArea).matNeuralRGR = cat_matNeuralRGR(:, catAreaID == idArea);
%     
% end

%% Merge only face patches
setAreaFP = 1:4;

matTS_FP.catChanID = cat(1, matTS_area(setAreaFP).setChanID); % catChanID;
matTS_FP.catSubjID = cat(1, matTS_area(setAreaFP).catSubjID); % catSubjID;
matTS_FP.setArea = {matTS_area(setAreaFP).nameArea}; % catSubjName;
matTS_FP.catAreaID = floor(matTS_FP.catSubjID./10); % catAreaID;
matTS_FP.matRaster = cat(2, matTS_area(setAreaFP).matRaster);
% matTS_FP.matSDF = cat(2, matTS_area(setAreaFP).matSDF);


% save the cocatenated time courses in different resolution
save(fullfile(directory.dataHome, sprintf('matRaster_Movie%s_allCells.mat', num2str(setMovie, '%d%d%d'))), ...
    'infoTS*', 'matTS*', 'paramTS', '-v7.3')



