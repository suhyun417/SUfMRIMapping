% renameFiles_Davida.m
%
% 2020/06/04 SHP
%        Updated to incorporate newly collected data in July 2018 in addition to the data from May 2018. 
%        Now the new excel file of the cell list that contains one
%        integrated list of cells across May and July dataset
%        (Davida_cell_list*_merged_SHP.xls). This excel contains the cell
%        identity across days (by Kenji/Elena) and evaluation of cell
%        quality (by SHP, based on reliability of responses across days to
%        fingerprinting stimuli and movies: 1 is good, 2 is so-so, 0 is
%        bad, -1 is strange (e.g. file is missing)
%
% 2018/08/17 SHP
%         Rename (and store) mat files containing Davida movie data collected by Elena/Kenji in May 2018, so that
%         they can be matched to other files previously collected from other monkeys
%         This code uses cell index in the excel file of
%         the cell list (e.g. Spice_cell_list_180120_SHP.xls), following the new version of "renameFiles_Spice.m"
%         This excel file contains the cell identity across days (by Kenji/Elena) and evaluation of
%         cell quality (by SHP, based on reliability of responses across days to
%         finger printing stimuli and movies: 1 is good, 2 is so-so, 0 is bad, -1 is strange (e.g. file is missing))

clear all;

nameSubjNeural = 'Dav'; 
directory.dataHome = '/procdata/parksh/_macaque';

% Directory info
directory.dataNeural = fullfile(directory.dataHome, nameSubjNeural); 
directory.source{1} = fullfile(directory.dataNeural, '_orgData', 'Davida180515_18_movie');
directory.source{2} = fullfile(directory.dataNeural, '_orgData', 'Davida180723-25_movie');
directory.destination = fullfile(directory.dataNeural);

% read the excel file of cell list
filename_xls =  fullfile(directory.dataNeural, '_orgData/Davida_cell_list_180515_180723_merged_SHP_movie123.xls'); %Davida_cell_list_180515_SHP.xls');
T = readtable(filename_xls, 'ReadRowNames', true);
setCol{1} = [1 2]; % set of column indices for the first (2018 May) dataset
setCol{2} = [3 4 5]; % set of column indices for the first (2018 May) dataset


for iS = 1:2
d_n =[]; clear listMatSUFile listSU_all listSUchannelID locB locC
d_n = dir(fullfile(directory.source{iS}, '*sig*.mat'));
[listMatSUFile{1:length(d_n), 1}] = deal(d_n.name); % Get the list of single unit files

diary(fullfile(directory.destination, sprintf('renameFiles_Davida_%s', date)))
diary on
% Copy files with new names
for iFile = 1:length(listMatSUFile)
    
    % find the movie id
    movID = char(regexp(listMatSUFile{iFile}, '\d*(?=sig)', 'match'));
    
    if ~ismember(str2double(movID), [1 2 3]) % only renaming data for movie 123
        continue;
    end

    curSU = strsplit(listMatSUFile{iFile}, {'sig', '.mat'}); %, 'CollapseDelimiters', true); 
    tempC = strsplit(curSU{2}, {'_', '-'}); %, 'CollapseDelimiters', true);
    
    indCell = [];
    for iCol = setCol{iS}
        indCell = find(strcmp(T{:, iCol}, strjoin(tempC(1:2), '_'))>0, 1);
        if ~isempty(indCell)
            break;
        end
    end
    
%     if sum(strcmp(tempC, 'x'))>0 % if any of the name contains "x" that means the unit was not sorted on a particular day
%         iCol = find(strcmp(tempC(2:3), 'x')==0); % go to the date that had the unit
%         tempLoc = find(strcmp(T.(iCol), strjoin(tempC([1, iCol+1]), '_'))>0);
%         indCell = intersect(find(ismissing(T(:, find(strcmp(tempC(2:3), 'x'))))), tempLoc);
%     else
%         indCell = find(strcmp(T.(1), strjoin(tempC(1:2), '_'))>0);
%         if length(indCell)>1
%             indCell = intersect(find(strcmp(T.(1), strjoin(tempC(1:2), '_'))>0), find(strcmp(T.(2), strjoin(tempC([1, 3]), '_'))>0));
%         end
%     end
    
    if isempty(indCell)
        continue;
    else
        if T.EVAL(indCell) > 0 % if the evaluation is okay
            newFileName = strcat(lower(nameSubjNeural), 'mov', movID, 'sig', sprintf('%02d', indCell), '.mat');
            source = fullfile(directory.source{iS}, listMatSUFile{iFile});
            destination = fullfile(directory.destination, newFileName);
            copyfile(source, destination)
            
            fprintf(1, '\n    %s is copied as %s', listMatSUFile{iFile}, destination);
        else % if the cell is not good, don't copy the file
            continue;
        end
    end
    
%     if length(tempC) < 3 % in case of the unit was clustered only from one day of recording session
%         indCell = find(strcmp(T.(iCol), strjoin(tempC(1:2), '_'))>0);
%         if isempty(indCell) 
%             indCell = find(strcmp(T.(iCol+1), strjoin(tempC(1:2), '_'))>0); % check next column
%         end
%     else
%         indCell = find(strcmp(T.(iCol), strjoin(tempC(1:2), '_'))>0); %find(cellfun(@isempty, strfind(T.(iCol), strjoin(tempC(1:2), '_')))<1); 
%     end
    
%     if T{indCell,5} > 0 % if the evaluation is okay
%         newFileName = strcat(lower(nameSubjNeural), 'mov', movID, 'sig', sprintf('%02d', indCell), '.mat');
%         source = fullfile(directory.source, listMatSUFile{iFile});
%         destination = fullfile(directory.destination, newFileName);
%         copyfile(source, destination)
%         
%         fprintf(1, '\n    %s is copied as %s', listMatSUFile{iFile}, destination);
%     else % if the cell is not good, don't copy the file
%         continue;
%     end    
end
end
diary off


% nRow = size(T, 1);
% for iRow = 1:nRow
%     if T.EVAL(iRow) < 1
%         continue;
%     end
%     
%     % determine which dataset this cell is from: column [1 2] is dataset 1,
%     % column [3 4 5] is dataset 2
%     if sum(find(~ismissing(T{iRow, 1:5}))<2)>0
%         iDataSet = 1;
%     else
%         iDataSet = 2;
%     end
% 
%     T{iRow, setCol{iDataSet}(1)}
%     
%     
% end
% 
% 


