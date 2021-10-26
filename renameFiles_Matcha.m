% renameFiles_Matcha.m
%
% 2017/10/23 SHP: modified from renameFiles_Spice.m
% Rename (and store) mat files containing Matcha movie data collected by Brian, so that
% they can be matched to other files previously collected from other
% monkeys for Park et al. (2017) Neuron paper
%
% 2018/08/29 SHP
% Select only the cells that have less than 30% of shared spikes with other
% cells. Cell selection details are in
% /procdata/parksh/Mat/Mat_SelectedCells_Movie.txt

dirDataNeural = '/procdata/parksh/Mat';
dirDataNeural_org = '/procdata/parksh/Mat/_orgData';

%%%%%%%%%%%%
% There needs to be some hard-coded values to fit to Russ dataset
% MatMov1to6_Movie%sData.mat (%s is 01 to 06)
%   has a struct with a same name: MatMov1to6_Movie01Data
%   MatMov1to6_Movie06Data(1) is for the first neuron
%       -- MatMov1to6_Movie06Data(1).chan = [67 1];
%       -- MatMov1to6_Movie06Data(1).mov_spikes
%           -- 1st col: movie number
%           -- 2nd col: spike times during movie presentation
%           -- 3rd col: msec level binning of the spikes during movie
%           -- 4th col: spike times during 60s pre-movie (baseline)
%           -- 5th col: msec level binning of the spikes during 60s pre-movie (baseline)
%           -- 6th col: spike times during 60s post-movie
%           -- 7th col: msec level binning of the spikes during 60s post-movie

setMovie = [1, 2, 3, 4, 5, 6];
numMovie = length(setMovie);

% Valid cells based on the shared-spike analysis: added 2018/08/29
% Selection details are in /procdata/parksh/Mat/Mat_SelectedCells_Movie.txt
indValidCell = [2     3     4     5     7     8     9    10    11    12    13    14    15    16    17    18 ...
    19    20    21    23    24    25    26    27    28    29    30];

abc = ['a':'z'];
for iMov = 1:numMovie
    nameFile = sprintf('MatMov1to6_Movie%02dData', setMovie(iMov));
    tempS = load(fullfile(dirDataNeural_org, nameFile));
    
    names = fieldnames(tempS);
    orgStruct = getfield(tempS, names{1});
    numCell = length(indValidCell); %size(orgStruct, 2);
    
    for iCell = 1:numCell
        
        indCell = indValidCell(iCell);
        
        % Renaming
        curCellID = sprintf('%03d%s', orgStruct(indCell).chan(1), abc(orgStruct(indCell).chan(2)));
        saveFileName = sprintf('mmov%dsig%s.mat', setMovie(iMov), curCellID);
        
        % Get the data
        dat.h.units = 'ms';
        dat.h.info.orgFileName = nameFile;
        
        indEmpty = cellfun('isempty', orgStruct(indCell).mov_spikes(:,2)); % get rid of empty arrays, if any
        dat.s = deal(orgStruct(indCell).mov_spikes(indEmpty<1,2));
        dat.s_prepost = deal(orgStruct(indCell).mov_spikes(indEmpty<1, [4, 6]));
        
        dat.h.info.blocks = orgStruct(indCell).blocks(indEmpty<1);
       
        % Save in new file
        save(fullfile(dirDataNeural, saveFileName), 'dat')
        fprintf(1, '\n Movie #%d, Cell #%d, Cell ID %s, saved in %s ...', setMovie(iMov), indCell, curCellID, saveFileName);
    end
end

