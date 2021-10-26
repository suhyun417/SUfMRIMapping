% renameFiles_Dexter.m
%
% 2018/10/26 SHP
% Rename (and store) mat files containing Dexter's Spice movie data collected by Aidan,
% so that they can be matched to other files previously collected from other
% monkeys for Park et al. (2017) Neuron paper

clear all;

dirDataNeural = '/procdata/parksh/Dex';
dirDataNeural_org = '/procdata/parksh/Dex/_orgData';
load(fullfile(dirDataNeural_org, 'FaceCellSpike.mat'))

% setMovie = [1, 2, 3, 4, 5, 6];
% numMovie = length(setMovie);

abc = ['a':'z'];
numCell = length(SpikeData);
% for iMov = 1:numMovie
%     nameFile = sprintf('MatMov1to6_Movie%02dData', setMovie(iMov));
%     tempS = load(fullfile(dirDataNeural_org, nameFile));
%
%     names = fieldnames(tempS);
%     orgStruct = getfield(tempS, names{1});
%     numCell = length(indValidCell); %size(orgStruct, 2);
for iCell = 1:numCell
    
    setMovie = SpikeData(iCell).MovieIDs;
    
    for iM = 1:length(setMovie)
        
        % Renaming
        curCellID = sprintf('%s%03d%s', SpikeData(iCell).info.Session, SpikeData(iCell).info.Channel, abc(SpikeData(iCell).info.Cell));
        saveFileName = sprintf('dexmov%dsig%s.mat', setMovie(iM), curCellID);
        
        % Get the data
        dat.h.units = 'ms';
        dat.h.info = SpikeData(iCell).info;
        dat.s = deal(SpikeData(iCell).SpikeTimes(:,iM));
        
        % Save in new file
        save(fullfile(dirDataNeural, saveFileName), 'dat')
        fprintf(1, '\n Cell #%d/%d, Cell ID %s, Movie #%d, saved in %s ...', iCell, numCell, curCellID, setMovie(iM), saveFileName);
    end
    
end

