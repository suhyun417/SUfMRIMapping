% saveMovieData_SU_SDF_eachArea.m
%
% 2019/05/22 SHP
% Modified from saveMovieData_SU_SDF_multipleSubjects.m
% Prepare time-series matrices in different resolution for each face patch
% across different animals

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

setSpiceDirectory = {'2016Nov_movie', '2018Jan_movie'};

setArea = {'AF', 'AM', 'AM+', 'ML'};
setNameSubjNeural{1} = {'Spi', 'Moc'}; % for AF, cells with fingerprinting results
setNameSubjNeural{2} = {'Mat', 'Was'}; % for AM
setNameSubjNeural{3} = {'Dan', 'Moc'}; % for AM+
setNameSubjNeural{4} = {'Dav'}; % Avalanche has movie 4-6, not 1-3

%% Prepare time-series matrices in different resolution for each area
for iArea = 1:2 % not yet for ML
    nameArea = setArea{iArea};
    
    setMovie = [1 2 3]; %[4 5 6]; %[1 2 3]; %[4 5 6]; %[1 2 3];
    nt = 375;
    typeSpiceData = 2; %1; % indicate whether it's 2016 data (1) or 2018 data (2)
    
    % MION function
    TR=2.4;
    k = gampdf([-40:TR:40],4,2);
    
    curSetNameSubjNeural = setNameSubjNeural{iArea};
    
    paramSDF_area.setMovie = setMovie;
    paramSDF_area.TR = TR;
    paramSDF_area.MIONfun = k;
    paramSDF_area.setNameSubjNeural = curSetNameSubjNeural;
    
    matSDF = struct([]);
    for iSubj = 1:length(curSetNameSubjNeural)
        
        % Load the data
        nameSubjNeural = curSetNameSubjNeural{iSubj};
        dirDataNeural = fullfile(directory.dataHome, nameSubjNeural);
        if sum(strcmpi(nameSubjNeural, {'spice', 'spi'})) && typeSpiceData == 1
            dirDataNeural = fullfile(dirDataNeural, '2016Nov_movie');
        end
        filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
        load(fullfile(dirDataNeural, filenameNeural))
        fprintf(1, '\nLoading single unit data of %s: %s ....\n ', nameSubjNeural, filenameNeural)
        
        [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
        validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
        
        if isempty(validC)
            continue;
        end
        
        if typeSpiceData == 1
            switch lower(nameSubjNeural)
                case 'spi'
                    excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
                    validC = setdiff(1:length(paramSDF.setCellIDs), excChanIndex)';
                otherwise
                    [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
                    validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
            end
        end
        
        if strcmpi(nameSubjNeural, 'moc') % Mochi has both AF and AM cells, so need to select cells from each area
            validC = find(~cellfun(@isempty, strfind(paramSDF.setCellIDs, nameArea)));
        end
        
        matSDF(iSubj).nameSubj = nameSubjNeural;
        matSDF(iSubj).area = nameArea; %setNameArea{iSubj};
        matSDF(iSubj).setCellIDs = paramSDF.setCellIDs(validC);
        matSDF(iSubj).indValidCell = validC;
        
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
        dirDataNeural_individualCellFile = dirDataNeural;
        if sum(strcmpi(nameSubjNeural, {'spice', 'spi'})) && typeSpiceData == 2
            dirDataNeural_individualCellFile = fullfile(dirDataNeural, setSpiceDirectory{typeSpiceData}); % '2018Jan_movie');
        end
        FR_dTfine = createCellRegressor_indMov_discreteTime(dirDataNeural_individualCellFile, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
            setMovie, 0.1); % in 10Hz (number of spikes for every 100ms)
        % concatenate across movies
        matFR_SU=[];
        for iUnit = 1:size(FR_dTfine,1)
            tempFR = cat(1, FR_dTfine(iUnit, :).mnFR);
            matFR_SU(:,iUnit) = tempFR;
        end
        matSDF(iSubj).matFR_SU = matFR_SU;
        
        % 4. Fine temporal resolution, normalized for each movie (z-scored)
        matFR_SU_norm=[];
        for iM = 1:length(setMovie)
            tempFR = cat(2, FR_dTfine(:, iM).mnFR);
            matFR_SU_norm = cat(1, matFR_SU_norm, zscore(tempFR)); % normalized time series for each movie and concatenate across movies
        end
        matSDF(iSubj).matFR_SU_norm = matFR_SU_norm;
        
    end
    
    % save the cocatenated time courses in different resolution
    save(fullfile(directory.dataHome, sprintf('matSDF_%s_Movie%s.mat', nameArea, num2str(setMovie, '%d%d%d'))), 'matSDF', 'paramSDF_area')
    
end