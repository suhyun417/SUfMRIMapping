% renameFiles_Mochi.m
%
% 2019/05/21 SHP
% Rename (and store) mat files containing Mochi movie data collected by Kenji & Elena W 
% in 1) October 2018 & 2) March 2019 (two different population), so that they can be matched to 
% other files previously collected from other monkeys.
%       - Modified from renameFiles_Dango.m
%       - Using cell index in the excel file of the cell list (e.g. Spice_cell_list_180120_SHP.xls)
%       - This excel file contains the cell identity across days (by Kenji/Elena) and evaluation of
%       cell quality (by SHP, based on reliability of responses across days to
%       finger printing stimuli and movies: 1 is good, 2 is so-so, 0 is bad, -1 is strange (e.g. file is missing))

clear all;

nameSubjNeural = 'Moc'; 
directory.dataHome = '/procdata/parksh/_macaque'; %'/procdata/parksh';

% Directory info
directory.dataNeural = fullfile(directory.dataHome, nameSubjNeural); 
% directory.source = fullfile(directory.dataNeural, '_orgData');
setSessionName = {'Mochi181023_movie', 'Mochi190313_movie'};
setXlsName = {'Mochi_cell_list_181023_SHP.xls', 'Mochi_cell_list_190313_SHP.xls'};
directory.destination = fullfile(directory.dataNeural);

% filename_xls{1} = '/procdata/parksh/_macaque/Moc/_orgData/Mochi_cell_list_181023_SHP.xls';
% filename_xls{2} = '/procdata/parksh/_macaque/Moc/_orgData/Mochi_cell_list_190313_SHP.xls'; 

% Area info for two sessions
setArea = {'AF', 'AM'};
setIndArea = zeros(128, 2); % channel x session
setIndArea(1:64, 1) = 2; % 181023: channel 1:64 is AM
setIndArea(65:128, 1) = 1; % 181023: channel 65:128 is AF
setIndArea(1:64, 2) = 1; % 190313: channel 1:64 is AF
setIndArea(65:128, 2) = 2; % 190313: channel 65:128 is AM

% Set of cells that need removal of the 7th trial from session 2
setUnitTrialRemoval = [4 10 12 14 17 18 20 23 29]; 

%% Copy and save with a new name: only for valid cells 
for iSession = 1:2
    directory.source = fullfile(directory.dataNeural, '_orgData', setSessionName{iSession});
    
    % read the excel file of cell list
    filename_xls =  fullfile(directory.dataNeural, '_orgData', setXlsName{iSession});
    T = readtable(filename_xls, 'ReadRowNames', true);
    
    d_n =[]; clear listMatSUFile listSU_all listSUchannelID locB locC
    d_n = dir(fullfile(directory.source, '*sig*.mat'));
    [listMatSUFile{1:length(d_n), 1}] = deal(d_n.name); % Get the list of single unit files
    
    diary(fullfile(directory.destination, sprintf('renameFiles_%s_%s', setSessionName{iSession}, date)))
    diary on
    % Copy files with new names
    for iFile = 1:length(listMatSUFile)
        
        % find the movie id
        movID = char(regexp(listMatSUFile{iFile}, '\d*(?=sig)', 'match'));
        %
        %     if ismember(str2double(movID), [4 5 6])
        %         iCol = 1; % look at the first column of the excel sheet
        %     elseif ismember(str2double(movID), [1 2 3])
        %         iCol = 3; % look at the third column of the excel sheet
        %     end
        
        curSU = strsplit(listMatSUFile{iFile}, {'sig', '.mat'}, 'CollapseDelimiters', true); % with Elena W's naming scheme with one underscore one dash
        tempC = strsplit(curSU{2}, {'_', '-'}, 'CollapseDelimiters', true);
%         iCol = 1; % by default, look at the first day
        if sum(strcmp(tempC, 'x'))>0 % if any of the name contains "x" that means the unit was not sorted on a particular day
            iCol = find(strcmp(tempC(2:3), 'x')==0); % go to the date that had the unit
            tempLoc = find(strcmp(T.(iCol), strjoin(tempC([1, iCol+1]), '_'))>0);
            indCell = intersect(find(ismissing(T(:, find(strcmp(tempC(2:3), 'x'))))), tempLoc);
        else
            indCell = find(strcmp(T.(1), strjoin(tempC(1:2), '_'))>0);
            if length(indCell)>1
                indCell = intersect(find(strcmp(T.(1), strjoin(tempC(1:2), '_'))>0), find(strcmp(T.(2), strjoin(tempC([1, 3]), '_'))>0));
            end
        end
        
        if T.EVAL(indCell) > 0 %T{indCell,5} > 0 % if the evaluation is okay
            nameArea = setArea{setIndArea(str2num(char(tempC(1))), iSession)};
            newFileName = strcat(lower(nameSubjNeural), 'mov', movID, 'sig', sprintf('%d%02d', iSession, indCell), nameArea, '.mat');
            source = fullfile(directory.source, listMatSUFile{iFile});
            destination = fullfile(directory.destination, newFileName);
            
             % check the struct name and copy
            tempS = load(source);    
            names = fieldnames(tempS);
            orgStruct = getfield(tempS, names{1});
            dat = orgStruct; clear orgStruct
            
            save(destination, 'dat')
            
            fprintf(1, '\n    %s is copied as %s', listMatSUFile{iFile}, destination);
            clear dat
            
%             if iSession==2 && ismember(indCell, setUnitTrialRemoval) %
%             decided not to do it because it has a marginal effect on the
%             final mean ts across trials
%                 load(destination, 'dat');
%                 % manually remove the 7th trial
%                 dat_org = dat;
%                 clear dat
%                 ttt = ones(1, length(dat_org.s));
%                 ttt(7) = 0;
%                 ttt = logical(ttt);
%                 dat.s = dat_org.s(ttt);
%                 dat.trial = dat_org.trial(ttt);
%                 dat.tank = dat_org.tank(ttt);
%                 dat.waveform = dat_org.waveform(ttt);
%                 dat.h = dat_org.h;
%                 dat.InfoRemovalTrial = '7th trial removed manually';
%                 clear dat_org
%                 save(destination, 'dat')
%                 fprintf(1, '\n                7th trial is manually removed from the cell #%02d and re-saved', indCell);
%             end
        else % if the cell is not good, don't copy the file
            continue;
        end
    end
    diary off
end
