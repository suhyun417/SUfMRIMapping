function [] = reconVoxCI(nameSubjNeural,nameSubjBOLD, iChan) %nameSubjBOLD, 
% Reconstruct swarm-generated results of confidence interval for each
% cell-voxel pair and save results for each cell
%
% 1. Load the confidence interval data for each set of 1000 voxels 
% 2. Merge in into new large struct "resultBS

% Set directories 
dirDataHome = '/data/parks20/procdata/NeuroMRI'; % Biowulf %'/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirCIFile = fullfile(dirDataNeural, 'corrMap_resultsCI_10hz');
% if ~exist(dirSaveFile, 'dir')
%     mkdir(dirSaveFile)
% end
% dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);


% % Directory for saving figures as graphic files
% dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';



%%
% Take care of possible string-double issue from compile
if ischar(iChan); eval(sprintf('iChan = %s;', iChan)); end;

setStartVox = 1:1000:81920;
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr') % get cell IDs from another file

% nChanSwarm = 12;
% setChan = idChan_start:1:idChan_start+nChanSwarm-1;

resultBS = struct();
% for iChan = 1:length(setChan)
    
%     idChan = setChan(iChan);
    
    resultBS.chanID = paramCorr.validChanID(iChan,:);
    
    countVox = 0;
    for ii = 1:length(setStartVox)
        
        idVox_start = setStartVox(ii);
        loadFileName = sprintf('corrcoeffCI_%s%s_cell%d_vox%d.mat', nameSubjNeural, nameSubjBOLD, iChan, idVox_start);
        
        load(fullfile(dirCIFile, loadFileName), 'VoxCI')
        
        resultBS.VoxCI(countVox+1:countVox+length(VoxCI)) = VoxCI;
        countVox = countVox + length(VoxCI);
    end
% end

saveFileName = sprintf('corrcoeffCI_%s%s_cell%d.mat',...
    nameSubjNeural, nameSubjBOLD, iChan);

save(fullfile(dirCIFile, saveFileName),'resultBS')
