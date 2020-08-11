% renameFigureFiles.m
%
% rename individual SU Map graphic files to match their cell ID

clear all;

nameSubjNeural = 'Spi';
nameSubjBOLD = 'Art';
dirPics = sprintf('/projects/parksh/NeuroMRI/_labNote/_figs/_CorrMap/SUmap_unmasked/%s/%s/', nameSubjBOLD, nameSubjNeural);

% load Cell IDs for renaming
% T = readtable(fullfile(dirPics, 'SpiceCells.txt'), 'Delimiter', '\t', 'ReadVariableNames', false, 'HeaderLines', 1, 'Format', '%d%s');
T = readtable(sprintf('/projects/parksh/NeuroMRI/analysis/%sCells.txt', nameSubjNeural), 'Delimiter', '\t', 'ReadVariableNames', false, 'HeaderLines', 1, 'Format', '%d%s');

% get the list of graphic files to rename
filestr = '*0p5f.0*.png'; % this is manually coded for each case
d = dir(fullfile(dirPics, filestr));

for iFile = 1:length(d)
    a = regexp(d(iFile).name, '\d*(?=.png)', 'match');
    if isempty(a)
        continue;
    end
    iCell =  find(T.(1)==str2double(a));
    cellID = T{iCell, 2};
    
    temp = strsplit(d(iFile).name, '.');
    fname_new = sprintf('%s_cell%s.png', temp{1}, char(cellID));
    
    movefile(fullfile(dirPics, d(iFile).name), fullfile(dirPics, fname_new));
    fprintf(1, '   :: Renaming file %s to %s :: \n', d(iFile).name, fname_new)
end
    
    