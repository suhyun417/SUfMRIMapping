

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


for iSubj = 1:length(setNameSubjNeural) %curSetNameSubjNeural)
    
    % Load the data
    nameSubjNeural = setNameSubjNeural{iSubj}; %curSetNameSubjNeural{iSubj};
    dirDataNeural = fullfile(directory.dataHome, nameSubjNeural);
    
    filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
    load(fullfile(dirDataNeural, filenameNeural))
    fprintf(1, '\nLoading single unit data of %s: %s ....\n ', nameSubjNeural, filenameNeural)
    
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
       
    % 10 millisecond binned time series for each trial
    sizeTimeBin_sec = 0.01;
    FR_10msBin = createCellRegressor_indMov_discreteTime(dirDataNeural, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
        setMovie, sizeTimeBin_sec); % in 100Hz (number of spikes for every 10ms)
    
end
    
    
    
    
    