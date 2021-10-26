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
% nameSubjNeural = 'Spi'; %'Tor';
% nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
% dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);


%% 
% setMovie = [1 2 3];
nMovie = length(setMovie);
for iM = 1:nMovie
    idMov = setMovie(iM);
    
    % pick a file: for each movie, all the cells contain the same eye data,
    % so pick the first cell
    W = what(dirDataNeural);
    k = strfind(W.mat, sprintf('mov%dsig', idMov));
    tempInd = cell2mat(cellfun(@isempty, k, 'uni', false))<1;
    