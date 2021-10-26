% renameFiles_Spice.m
%
% Rename (and store) mat files containing Spice movie data collected by Kenji, so that
% they can be matched to other files previously collected from other
% monkeys
%
% 2018/08/16 SHP
%       Modified to rename files using cell index in the excel file of
%       the cell list (e.g. Spice_cell_list_180120_SHP.xls). This excel file
%       contains the cell identity across days (by Kenji/Elena) and evaluation of
%       cell quality (by SHP, based on reliability of responses across days to
%       finger printing stimuli and movies: 1 is good, 2 is so-so, 0 is bad)


clear all;

nameSubjNeural = 'Spi';
directory.dataHome = '/procdata/parksh';

% Directory info
directory.dataNeural = fullfile(directory.dataHome, nameSubjNeural); % '/procdata/parksh/Spi';
% dirDataNeural_org = fullfile(dirDataNeural, '_orgData'); % '/procdata/parksh/Spi/_orgData';
directory.source = fullfile(directory.dataNeural, '_orgData', 'Spice180120_26_movie');
directory.destination = fullfile(directory.dataNeural,  '2018Jan_movie');

% nameSession_source = 'Spice180120_26_movie';
% nameSession_destination = '2018Jan_movie';

% read the excel file of cell list
filename_xls =  fullfile(directory.dataNeural, '_orgData/Spice180120_other/Spice_cell_list_180120_SHP.xls');
T = readtable(filename_xls, 'ReadRowNames', true);


d_n =[]; clear listMatSUFile listSU_all listSUchannelID locB locC
d_n = dir(fullfile(directory.source, '*sig*.mat'));
[listMatSUFile{1:length(d_n), 1}] = deal(d_n.name); % Get the list of single unit files

diary(sprintf('renameFiles_Spice_%s', date))
diary on
% Copy files with new names
for iFile = 1:length(listMatSUFile)
    
    % find the movie id
    movID = char(regexp(listMatSUFile{iFile}, '\d*(?=sig)', 'match'));
    
    % find the matched cell index from the excel sheet
    if ismember(str2double(movID), [4 5 6])
        iCol = 1; % look at the first column of the excel sheet
    elseif ismember(str2double(movID), [1 2 3])
        iCol = 3; % look at the third column of the excel sheet
    end
    
    curSU = regexp(listMatSUFile{iFile}, '(?<=sig)\w*', 'match');
    tempC = strsplit(curSU{1}, '_');
    if length(tempC) < 3 % in case of the unit was clustered only from one day of recording session
        indCell = find(strcmp(T.(iCol), strjoin(tempC(1:2), '_'))>0);
        if isempty(indCell) 
            indCell = find(strcmp(T.(iCol+1), strjoin(tempC(1:2), '_'))>0); % check next column
        end
    else
        indCell = find(strcmp(T.(iCol), strjoin(tempC(1:2), '_'))>0); %find(cellfun(@isempty, strfind(T.(iCol), strjoin(tempC(1:2), '_')))<1); 
    end
    
    if T{indCell,5} > 0 % if the evaluation is okay
        newFileName = strcat(lower(nameSubjNeural), 'mov', movID, 'sig', sprintf('%02d', indCell), '.mat');
        source = fullfile(directory.source, listMatSUFile{iFile});
        destination = fullfile(directory.destination, newFileName);
        copyfile(source, destination)
        
        fprintf(1, '\n    %s is copied as %s', listMatSUFile{iFile}, destination);
    else % if the cell is not good, don't copy the file
        continue;
    end    
end
diary off



%     listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')';
%     [listSUchannelID, locB, locC] = unique(listSU_all);

% for iSession = 1:length(nameSession)
%     
%     channels=[];
%     abc = ['a':'z'];
%     for iCell = 1:length(listSUchannelID)
%         tempC = strsplit(listSUchannelID{iCell}, '_');
%         channels(iCell,1) = str2double(tempC{1});
%     end
%     
%     
%     listSUchannelID_newID=cell(size(listSUchannelID));
%     for iCell = 1:length(listSUchannelID)
%         tempC = strsplit(listSUchannelID{iCell}, '_');
%         charabc='a'; % default
%         if sum(channels==str2double(tempC{1}))>1 % if there are more than one unit from this channel
%             temploc = find(channels==str2double(tempC{1}));
%             orderUnit = find(temploc == iCell);
%             charabc=abc(orderUnit);
%         end
%         newID = sprintf('%03d%s', str2double(tempC{1}), charabc);
%         listSUchannelID_newID{iCell} = newID;
%     end
%     
%     % Copy files with new names
%     for iFile = 1:length(listMatSUFile)
%         if ~strcmp(listMatSUFile{iFile}, 'mov') % if the original name doesn't contain "mov" string
%             % find the movie id
%             movID = char(regexp(listMatSUFile{iFile}, '\d*(?=sig)', 'match'));
%             newFileName = strcat(lower(nameSubjNeural(1)), 'mov', movID, 'sig', listSUchannelID_newID{locC(iFile)}, '.mat');
%         else
%             newFileName = strrep(listMatSUFile{iFile}, listSU_all{iFile}, listSUchannelID_newID{locC(iFile)});
%         end
%         source = fullfile(dirDataNeural_org, nameSession{iSession}, listMatSUFile{iFile});
%         destination = fullfile(dirDataNeural,  '2018Jan_movie', newFileName);
%         copyfile(source, destination)
%     end
%     
% end