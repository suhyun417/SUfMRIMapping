% genFig_multipleFP_Revision_crossorrelationCells.m
%
% 2021/11/1 SHP
% Compute cross-correlations across neurons in each area using
% highpass-filtered time series

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

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was', 'Dav'};

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

setMovie = [1 2 3]; %[4 5 6]; %[1 2 3]; %[4 5 6]; %[1 2 3];
nt = 375;

% MION function
TR=2.4;
k = gampdf([-40:TR:40],4,2);

%     paramSDF.setNameSubjNeural = curSetNameSubjNeural;
% infoTS_subj = struct([]); matTS_subj = struct([]);
setArea = {'AF', 'AM', 'AAM', 'ML', 'NFP'};
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
       
    indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
    % 1. No MION, TR resolution (375 time points)
    matFR_TR = NaN(nt, length(validC));
    for iChan = 1:length(validC)
        tempFR = [];
        tempFR = cat(1, S(validC(iChan), indMovieNeuron).mnFR);
        matFR_TR(:,iChan) = tempFR;
    end
    matTS_subj(iSubj).matFR_TR = matFR_TR;
    
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
    matTS_subj(iSubj).matNeuralRGR = matNeuralRGR;
    
    % 3. Fine temporal resolution, normalized (z-scored) for each movie
    sizeTimeBin_sec = 0.1;
    FR_dTfine_10hz = createCellRegressor_indMov_discreteTime(dirDataNeural, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
        setMovie, sizeTimeBin_sec); % in 10Hz (number of spikes for every 100ms)
    % concatenate across movies
    matFR_SU=[];
    for iUnit = 1:size(FR_dTfine_10hz,1)
        tempFR = cat(1, FR_dTfine_10hz(iUnit, :).mnFR);
        matFR_SU(:,iUnit) = tempFR;
    end
    matTS_subj(iSubj).matFR_SU_10hz = matFR_SU;