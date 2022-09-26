function [indDataMovie, listSUchannelID, listMovie] = genDataMatrix_SU(nameSubj, flagFig)
% Generate data matrix
%
% Make a data matrix that contains logical values for data existence of a
% particular cell (channel) for a movie

ss = pwd;
dirData = '/nifvault/procdata/parksh/_macaque';
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirData = '/Volumes/NIFVAULT/procdata/parksh/_macaque';   
end
dirDataNeural = fullfile(dirData, nameSubj); 
% if sum(strcmpi(nameSubj, {'spice', 'spi'}))
%     dirDataNeural = fullfile(dirData, nameSubj, '2018Jan_movie');
% end
d_n = dir(fullfile(dirDataNeural, '*sig*.mat')); %'*LFP*')); %'*sig*.mat'));
[listSUFile{1:length(d_n), 1}] = deal(d_n.name);

% list of channel IDs and movie IDs of each file
listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')'; % list of channels
listSUchannelID = unique(listSU_all); % only unique channels
listMovie = sort(str2num(char(unique(regexp(cat(2,d_n.name), '\d*(?=sig)', 'match')')))); % list of movies

% Since not every channel has data for every movie, 
% the "data presence matrix (Channels x MovieIDs)" is generated here
% for which channel has which movie's data: 1 for data, 0 for no data

% listMovie = unique(listMov_all);
indDataMovie=[];
for iCh=1:length(listSUchannelID) % go throucgh channel-by-channel     
     tempListMov{iCh} = regexp([listSUFile{strcmp(listSUchannelID{iCh}, listSU_all)}], '\d*(?=sig)', 'match');
     % Get indices for common movies across cells
     indDataMovie(iCh, :) = ismember(listMovie, str2num(char(tempListMov{iCh}))'); % ismember(listMovie, tempListMov{iCh}); %ismember(listMovie, str2num(char(tempListMov{iCh}))'); % data presence matrix (Channels x MovieIDs) 
end

if flagFig
    figure;
    imagesc(indDataMovie*(-1))
    colormap(gray)
end

end