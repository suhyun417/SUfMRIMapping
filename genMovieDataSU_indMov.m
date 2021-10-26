function [] = genMovieDataSU_indMov() %, saveFileName_matSDF, saveFileName)

%
% Save single unit data for individual movies for a particular subject under the individual data folder
%
% 11/11/14, SHP

% Subject name and directory
% nameSubjNeural = 'Tor'; %'Rho';
% nameSubjBOLD = 'Art';
% dirDataNeural = [dirData, nameSubjNeural, '/'];
fprintf(1, '\nSave single unit time series (spike density function) for individual movies separately\n');
nameSubjNeural = input('\nName of subject? (e.g. Tor, Rho):', 's');

if sum(strcmpi(nameSubjNeural, {'toroid', 'tor', '', ''}))
    nameSubjNeural = 'Tor';
elseif sum(strcmpi(nameSubjNeural, {'rhombus', 'rho', ''}))
    nameSubjNeural = 'Rho';
elseif sum(strcmpi(nameSubjNeural, {'Ziggy', 'Sieglinde', 'Si', 'Sieggie', 'Sig', 'AZ26', 'AZ2'}))
    nameSubjNeural = 'Sig';
end

dirData = '/procdata/parksh/';
dirDataNeural = [dirData, nameSubjNeural, '/'];
saveFileName = fullfile(dirDataNeural, sprintf('%s_movieTS_SU_indMov.mat', nameSubjNeural));

fprintf(1, '\n Save single unit time series for subject "%s"', nameSubjNeural);
fprintf(1, '\n Data will be saved as "%s" \n ', saveFileName);
fprintf(1, '\nPress enter to proceed \n'); input('')


%%
% Get the movie list from single unit data
d_n = dir(fullfile(dirDataNeural, '*sig*.mat'));
[listMatSUFile{1:length(d_n), 1}] = deal(d_n.name); % Get the list of single unit files 
% listMatSUFile = regexp(cat(2, d_n.name), '\w*.mat', 'match')'; %regexp(cat(2, d_n.name), '(?<=mat)\w*', 'match')';
% listSUchannelID = unique(regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match'));
listMovSU = unique(regexp(cat(2,d_n.name), '\d*(?=sig)', 'match')); % Movie IDs of the single units in cellstr

% % Select common movie IDs between single units and BOLD (or define the movies)
% setMovIDs = input('Movie IDs (e.g. [1 2 3])? (default: common movies in single units and fMRI)');
% if isempty(setMovIDs)
%     setMovIDs = intersect(listMovBOLD, str2num(char(listMovSU))');
% end
% % setMovIDs = [13 14 15]; %intersect(listMovBOLD, str2num(char(listMovSU))');
% [~, locMovBOLD] = intersect(listMovBOLD, setMovIDs);
% % [setMovIDs, locMovBOLD] = intersect(listMovBOLD, str2num(char(listMovSU))'); 

% Get the list of single units & their movie
listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')';
listSUchannelID = unique(listSU_all);

% indDataMov=[];
% for iCh=1:length(listSUchannelID) % go throucgh channel-by-channel     
%      tempListMov{iCh} = regexp([listMatSUFile{strcmp(listSUchannelID{iCh}, listSU_all)}], '\d*(?=sig)', 'match');
%      % Get indices for common movies across cells
%      indDataMov(iCh, :) = ismember(setMovIDs, str2num(char(tempListMov{iCh}))'); % binary matrix for data presence (Channels x MovieIDs) 
% end

% setCellIDs = listSUchannelID(sum(indMov,2)==length(setMovIDs)); % For now, let's do it with cells that have data for every movies



%% Get mean SDF and regressors for each single unit
setCellIDs = listSUchannelID; % should be cell array of string
setMovIDs = sort(str2num(char(listMovSU))');
sizeTimeBin_sec = 2.4; %TR
S = createCellRegressor_indMov_discreteTime(dirDataNeural, setCellIDs, setMovIDs, sizeTimeBin_sec);
% -- S(iChan, iMov): for each cell & each movie (S.cellID), 
%       1) matFR: for each movie (cell), matSDF for each trial, 2) mnFR: mean SDF across trials, 

paramSDF.sizeTimeBin_sec = sizeTimeBin_sec;

% S = createCellRegressor_indMov(dirDataNeural, setCellIDs, setMovIDs, indDataMov); %, flagConcat); 
% % -- S(iChan, iMov): for each cell & each movie (S.cellID), 
% %       1) matsdf: for each movie (cell), matSDF for each trial, 2) mnsdf: mean SDF across trials, 
% %       3) catmnsdf: concatenated mean SDF across movies
% %       4) rgrsMION: MION function convolved regressors (before resampling)
% %       5) rgrsMION_resample: MION function convolved regressors after
% %       resampling (in accordance with the number of frames in fMRI data)


%% save File
save(saveFileName, 'S', 'paramSDF')
fprintf(1, 'Time series are saved in %s\n', saveFileName)
