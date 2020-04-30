function [] = genMovieDataBOLD_indMov() %, saveFileName_matSDF, saveFileName) %[dataBOLD] = genMovieDataBOLD_indMov(nameSubjBOLD) %,

%
% Save fMRI data for individual movies for a particular subject under the individual data folder
%
% 11/11/14, SHP

% Subject name and directory
% nameSubjNeural = 'Tor'; %'Rho';
% nameSubjBOLD = 'Art';
% dirDataNeural = [dirData, nameSubjNeural, '/'];
fprintf(1, '\nSave fMRI time series for individual movies separately\n');
nameSubjBOLD = input('\nName of subject? (e.g. Art, Ava):', 's');

dirData = '/procdata/parksh/';
dirDataBOLD = [dirData, nameSubjBOLD, '/'];
saveFileName = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

% Get the movie list from the fMRI data
fprintf(1, '\Loading BOLD time series...\n');

d_b = dir(fullfile(dirDataBOLD, '*BOLD.mat')); %W = what(dirDataBOLD);
load(fullfile(dirDataBOLD, d_b(1).name)); % for now it's only one file
setMovIDs = MCD.unimov; % Movie IDs of the current BOLD


%% First, split fMRI responses into each individual movie

fprintf(1, 'Number of movies in original time series: %d\n', length(setMovIDs))

dataBOLD = MCD; % Get necessary params for further DSP in BlockAna: for now, just all the params
dataBOLD = rmfield(dataBOLD, 'catmvoltc'); %{'catmvoltc', 'catmdgz', 'catimgdat'}); % remove the previous variable of concatenated time series

lengthIndMovie_TR = dataBOLD.max_tr;
% lengthIndMovie_sec = lengthIndMovie_TR*dataBOLD.TR; % 300;

% Get the necessary part of time courses, corresponding to our (common) movie set
% if ~issorted(setMovIDs) % just to be sure
%     setMovIDs = sort(setMovIDs);
% end

% countTR = 0; tc=[];
% Split concatenated time series to individual cells for individual movies
for iMov = 1:length(setMovIDs) 
    fprintf(1, 'Movie ID: %d (%d/%d), Location in original BOLD: %d\n', setMovIDs(iMov), iMov, length(setMovIDs), setMovIDs(iMov));
%     tc(:,:,:,countTR+1:countTR+lengthIndMovie_TR) = MCD.catmvoltc(:,:,:,lengthIndMovie_TR*(locMovBOLD(iMov)-1)+1:lengthIndMovie_TR*locMovBOLD(iMov)); % should be optimized for different set of movies
    dataBOLD.mvoltc{iMov} = MCD.catmvoltc(:,:,:,lengthIndMovie_TR*(setMovIDs(iMov)-1)+1:lengthIndMovie_TR*setMovIDs(iMov)); % should be optimized for different set of movies
    %     countTR = size(tc, 4);
end
% dataBOLD.catmvoltc = tc;
% fprintf(1, 'Time course for %d movies has been saved in dataBOLD.tc\n', length(setMovIDs));
% fprintf(1, 'The number of TRs of dataBOLD.catmvoltc should be %d and is %d\n', ...
%     lengthIndMovie_TR*length(setMovIDs), size(dataBOLD.catmvoltc, 4));

dataBOLD.unimov = setMovIDs;

clear MCD %tc %SigData

% save File
save(fullfile(dirDataBOLD, saveFileName), 'dataBOLD')
fprintf(1, 'Time series is saved in %s\n', fullfile(dirDataBOLD, saveFileName))


% % Select common movie IDs between single units and BOLD (or define the movies)
% setMovIDs = input('Movie IDs (e.g. [1 2 3])? (default: common movies in single units and fMRI)');
% if isempty(setMovIDs)
%     setMovIDs = intersect(listMovBOLD, str2num(char(listMovSU))');
% end
% % setMovIDs = [13 14 15]; %intersect(listMovBOLD, str2num(char(listMovSU))');
% [~, locMovBOLD] = intersect(listMovBOLD, setMovIDs);
% % [setMovIDs, locMovBOLD] = intersect(listMovBOLD, str2num(char(listMovSU))'); 
% 
% % flagConcat = input('Do you want to concatenate timeseries across movies? (1:yes, 0:no):');
% 
% % Get the list of single units & their movie
% listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')';
% listSUchannelID = unique(listSU_all);
% 
% indDataMov=[];
% for iCh=1:length(listSUchannelID) % go throucgh channel-by-channel     
%      tempListMov{iCh} = regexp([listMatSUFile{strcmp(listSUchannelID{iCh}, listSU_all)}], '\d*(?=sig)', 'match');
%      % Get indices for common movies across cells
%      indDataMov(iCh, :) = ismember(setMovIDs, str2num(char(tempListMov{iCh}))'); % data presence matrix (Channels x MovieIDs) 
% end
% 
% % % For now, let's do it with cells that have data for every movies
% % setCellIDs = listSUchannelID(sum(indMov,2)==length(setMovIDs));
% setCellIDs = listSUchannelID;
% 
% % setCellIDs = {'006a'; '009a'; '010a'; '012a'; '013a'}; % Rhombus %{'003a';'005a';'007a';'009a';'013b';'014a'}; % Sig: '013a' has very low spikes
% % setMovIDs  = [7 8 9 10 11 12]; % Rhom %[7 8 9]; % Sig % [1 2 3 7 8 9];





