% saveMovieDataSU_indMov_10fps.m
%
% 2018/02/28 SHP
% Save single unit time courses (spike counts within unoverlapping time window) for each movie in 10fps
% To apply Xiaomin's kernels to our neural data

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
    
% % Add necessary toolbox 
% addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution

% Set directories 
setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi'};
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
% dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);
 

% % load the masks (movie-driven & brain-only mask)
% load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp'); %, 'brainMask_BlockAna3D');


%% Load the corr map and select valid channel: subject by subject
numSubject = size(setNameSubjNeural, 2);
% matR_SU_all = [];
matTS{1} = [];
matTS{2} = [];
matTS{3} = [];

for iSubj = 1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123_new.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr');
    
    % 2017/01/19: Excluded some units, in case of Spice's data which hasn't been went through this procedure yet
    switch lower(nameSubjNeural)
        case 'spi'
            excChanIndex = [10 13 22 27 30 49]; % cells were not same acrossd two days
            validChanIndex = setdiff(paramCorr.validChanIndex, excChanIndex);
            validChanID = cat(2, paramCorr.validChanID(validChanIndex,:), num2str(ones(size(validChanIndex)).*iSubj)); %paramCorr.validChanID(validChanIndex_clustering,:);
        otherwise
            validChanIndex = (1:length(paramCorr.validChanIndex))';
            validChanID = cat(2, paramCorr.validChanID(validChanIndex,:), num2str(ones(size(paramCorr.validChanIndex)).*iSubj)); %paramCorr.validChanID;
    end
    
    paramSDF(iSubj).nameSubj = nameSubjNeural;
    paramSDF(iSubj).validChanIndex = validChanIndex;
    paramSDF(iSubj).validChanID = validChanID;   
    
    
    setCellIDs = cellstr(paramCorr.validChanID(validChanIndex,:)); %listSUchannelID; % should be cell array of string
    setMovIDs = [1 2 3]; %sort(str2num(char(listMovSU))');
    sizeTimeBin_sec = 0.1; %0.033; %2.4; %TR
    
    for iMov = 1:3
        clear S 
        S = createCellRegressor_indMov_discreteTime(dirDataNeural, setCellIDs, iMov, sizeTimeBin_sec);
        matTS{iMov} = cat(2, matTS{iMov}, cat(2, S.mnFR));
    end
    
    paramSDF(iSubj).setMovIDs = setMovIDs;
    paramSDF(iSubj).sizeTimeBin_sec = sizeTimeBin_sec;
   
end

%% save File
saveFileName = fullfile(dirDataNeural, sprintf('%s_movieTS_SU_movie123_10fps.mat',  cell2mat(setNameSubjNeural)));
save(saveFileName, 'matTS', 'paramSDF')




