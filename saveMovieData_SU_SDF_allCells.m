% saveMovieData_SU_SDF_allCells.m
%
% 2019/06/20 SHP
% Modified from saveMovieData_SU_SDF_eachArea.m
% Prepare time-series matrices in different resolution for each face patch
% across different animals for all the neurons

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
        directory.projects = '/projects';
        directory.procdata = '/procdata';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/library';
        addpath(fullfile(directory.library, 'matlab_utils'));
    end
end

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was'};
% setNameArea = {'AF', 'AF', 'AF', 'AF', 

% setArea = {'AF', 'AM', 'AM+', 'ML'};
% setNameSubjNeural{1} = {'Tor', 'Rho', 'Sig', 'Spi', 'Moc'}; % for AF: 11 12 13 14 15
% setNameSubjNeural{2} = {'Mat', 'Was'}; % for AM: 21 22
% setNameSubjNeural{3} = {'Dan', 'Moc'}; % for AM+: 31 32
% setNameSubjNeural{4} = {'Dav'}; % Avalanche has movie 4-6, not 1-3

%
% nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
% dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
% load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'matR_SU', 'paramCorr');

%% Prepare time-series matrices in different resolution for all cells
% for iArea = 1:3 % not yet for ML
%     nameArea = setArea{iArea};

setMovie = [1 2 3]; %[4 5 6]; %[1 2 3]; %[4 5 6]; %[1 2 3];
nt = 375;
%     typeSpiceData = 2; %1; % indicate whether it's 2016 data (1) or 2018 data (2)

% MION function
TR=2.4;
k = gampdf([-40:TR:40],4,2);

% curSetNameSubjNeural = setNameSubjNeural{iArea};
% setSubjIDNumber = 10*iArea + [1:length(curSetNameSubjNeural)];

%     paramSDF.setNameSubjNeural = curSetNameSubjNeural;

matSDF = struct([]);
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
        case 'dan' %AM+
            subjID = 31;
        case 'moc' % AF and AM+
            subjID = [15 32];
        case 'was' %AM
            subjID = 22;
    end
    
    [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
    validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
    
    if isempty(validC)
        continue;
    end   
    
    matSDF(iSubj).nameSubj = nameSubjNeural;
%     matSDF(iSubj).area = setNameArea{iSubj};
    matSDF(iSubj).setCellIDs = paramSDF.setCellIDs(validC);
    matSDF(iSubj).indValidCell = validC;
    matSDF(iSubj).setCellIDs_subjID = ones(size(matSDF(iSubj).setCellIDs)).*subjID(1);
    if strcmpi(nameSubjNeural, 'moc') % Mochi has both AF and AM cells, so need to select cells from each area
        locAM = ~cellfun(@isempty, strfind(paramSDF.setCellIDs, 'AM'));
        matSDF(iSubj).setCellIDs_subjID(locAM) = subjID(2);
    end
    
    indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
    % 1. No MION, TR resolution (375 time points)
    matFR_TR = NaN(nt, length(validC));
    for iChan = 1:length(validC)
        tempFR = [];
        tempFR = cat(1, S(validC(iChan), indMovieNeuron).mnFR);
        matFR_TR(:,iChan) = tempFR;
    end
    matSDF(iSubj).matFR_TR = matFR_TR;
    
    % 2. MION convolved, TR resolution
    matNeuralRGR = NaN(nt, length(validC));
    for iChan = 1:length(validC)
        neuralrgrs=[];
        for iMov = 1:length(indMovieNeuron)
            curNeuralTC = S(validC(iChan), indMovieNeuron(iMov)).mnFR(8:125); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
            curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
            curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
            curNeuralTC = cat(2, NaN(1,7), curNeuralTC); %curNeuralTC(1:7) = NaN;
            neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
        end
        matNeuralRGR(:,iChan) = neuralrgrs';
    end
    matSDF(iSubj).matNeuralRGR = matNeuralRGR;
    
    % 3. Fine temporal resolution
%     dirDataNeural_individualCellFile = dirDataNeural;
%     if sum(strcmpi(nameSubjNeural, {'spice', 'spi'})) && typeSpiceData == 2
%         dirDataNeural_individualCellFile = fullfile(dirDataNeural, setSpiceDirectory{typeSpiceData}); % '2018Jan_movie');
%     end
    FR_dTfine = createCellRegressor_indMov_discreteTime(dirDataNeural, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
        setMovie, 0.1); % in 10Hz (number of spikes for every 100ms)
    % concatenate across movies
    matFR_SU=[];
    for iUnit = 1:size(FR_dTfine,1)
        tempFR = cat(1, FR_dTfine(iUnit, :).mnFR);
        matFR_SU(:,iUnit) = tempFR;
    end
    matSDF(iSubj).matFR_SU_10hz = matFR_SU;
    
    % 4. Fine temporal resolution, normalized for each movie (z-scored)
    matFR_SU_norm=[];
    for iM = 1:length(setMovie)
        tempFR = cat(2, FR_dTfine(:, iM).mnFR);
        matFR_SU_norm = cat(1, matFR_SU_norm, zscore(tempFR)); % normalized time series for each movie and concatenate across movies
    end
    matSDF(iSubj).matFR_SU_10hz_norm = matFR_SU_norm;
    
end

% end

paramSDF.setMovie = setMovie;
paramSDF.TR = TR;
paramSDF.MIONfun = k;
paramSDF.setSubjID = 'for AF: 11 12 13 14 15; for AM: 21 22; for AM+: 31 32';

% save the cocatenated time courses in different resolution
save(fullfile(directory.dataHome, sprintf('matSDF_Movie%s_allCells.mat', num2str(setMovie, '%d%d%d'))), 'matSDF', 'paramSDF')
