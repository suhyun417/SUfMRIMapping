% renameFiles_Dango_noFacePatch.m
%
% Rename (and store) mat files containing Dango movie data collected by Kenji in January 2018, 
% esp. the right hemisphere outside face patch population, so that
% they can be matched to other files previously collected from other monkeys
% modified from renameFiles_Dango.m
%
% 2020/06/25 SHP
%       Modified to rename files using cell index in the excel file of
%       the cell list (Dango_cell_list_180123_2_SHP.xls)
%       This excel file contains the cell identity across days (by Kenji/Elena) and evaluation of
%       cell quality (by SHP, based on reliability of responses across days to
%       finger printing stimuli and movies: 1 is good, 2 is so-so, 0 is bad, -1 is strange (e.g. file is missing))

clear all;

nameSubjNeural = 'Dan'; 
directory.dataHome = '/procdata/parksh/_macaque';

% Directory info
directory.dataNeural = fullfile(directory.dataHome, nameSubjNeural); 
directory.source = fullfile(directory.dataNeural, '_orgData', 'Dango180123_26_movie_noFacePatch');
directory.destination = fullfile(directory.dataNeural);

% read the excel file of cell list
filename_xls =  fullfile(directory.dataNeural, '_orgData/Dango_cell_list_180123_2_SHP.xls');
T = readtable(filename_xls, 'ReadRowNames', true);


d_n =[]; clear listMatSUFile listSU_all listSUchannelID locB locC
d_n = dir(fullfile(directory.source, '*sig*.mat'));
[listMatSUFile{1:length(d_n), 1}] = deal(d_n.name); % Get the list of single unit files

diary(fullfile(directory.destination, sprintf('renameFiles_Dango_noFacePatch_%s', date)))
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
        newFileName = strcat(lower(nameSubjNeural), 'mov', movID, 'sig', sprintf('NFP%02d', indCell), '.mat');
        source = fullfile(directory.source, listMatSUFile{iFile});
        destination = fullfile(directory.destination, newFileName);
        copyfile(source, destination)
        
        fprintf(1, '\n    %s is copied as %s', listMatSUFile{iFile}, destination);
    else % if the cell is not good, don't copy the file
        continue;
    end    
end
diary off




% nameSubjNeural = 'Dan';
% dirDataHome = '/procdata/parksh';
% 
% dirDataNeural = fullfile(dirDataHome, nameSubjNeural); % '/procdata/parksh/Spi';
% dirDataNeural_org = fullfile(dirDataNeural, '_orgData'); % '/procdata/parksh/Spi/_orgData';
% 
% nameSession{1} = 'Dango180123_26_movie';
% 
% for iSession = 1:length(nameSession)
%     
%     d_n = dir(fullfile(dirDataNeural_org, nameSession{iSession}, '*sig*.mat'));
%     [listMatSUFile{1:length(d_n), 1}] = deal(d_n.name); % Get the list of single unit files
%     listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')';
%     [listSUchannelID, locB, locC] = unique(listSU_all);
%     
%     % Renaming
%     channels=[];
%     abc = ['a':'z'];
%     for iCell = 1:length(listSUchannelID)
%         tempC = strsplit(listSUchannelID{iCell}, '_');
%         channels(iCell,1) = str2double(tempC{1});
%     end
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
%         
%         % excluding some trials that don't have spiking data
%         load(fullfile(dirDataNeural_org, nameSession{iSession}, listMatSUFile{iFile}))
%         locInvalidTrial = cellfun(@isempty, dat.s); 
%         if sum(locInvalidTrial) > 0 % if there's any trial without spikes
%             sprintf('Filename: %s contains %d invalid trials \n', listMatSUFile{iFile}, sum(locInvalidTrial))
%             dat.locInvalidTrial = locInvalidTrial;
%             dat.s_org = dat.s; % keep the original spike timing
%             s_new = deal(dat.s(logical(-(locInvalidTrial)+1)));
%             dat.s = s_new; % exclude the invalid trials with no spiking
%             clear s_new
%         end
%         source = fullfile(dirDataNeural_org, nameSession{iSession}, listMatSUFile{iFile});
%         destination = fullfile(dirDataNeural, newFileName);
%         copyfile(source, destination)
%     end
%     
% end