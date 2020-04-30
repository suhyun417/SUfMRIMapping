% getSplitHalfCorr.m
%
% Compute split-half correlation between two sets of trials of the same
% neuron and/or of the different neurons
% For now, this code follows the method used in McMahon et al., 2015 JN
% to compare the results across McMahon et al's results
% 2017/11/27 SHP

clear all;

dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

nameSubjNeural = 'Tor'; %'Mat'; % 'Spi'; %'Dan'; %'Dav'; %'Dan'; %'Spi'; %'Mat';
% Load data files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);

filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
load(fullfile(dirDataNeural, filenameNeural), 'paramSDF') % to get the cell ID

setCellID = paramSDF.setCellIDs;
setMovie = [1 2 3]; %[1 2 3 4 5 6]; %[1 2 3];
% sampleCell = {'065a', '075a', '118a', '082a'}; %'/procdata/parksh/Tor/
% win = [0 300]; % seconds
paramSplitHalfCorr.setMovie = setMovie;

% Compute the SDF 
paramSplitHalfCorr.sigma = 100; % sigma (for SDF computation) in msec
if sum(strcmpi(nameSubjNeural, {'spice', 'spi'}))
    dirDataNeural = [dirDataHome, nameSubjNeural, '/', '2018Jan_movie'];
end
S = createCellRegressor_indMov(dirDataNeural, setCellID, setMovie, paramSplitHalfCorr.sigma); %{'065a', '082a', '089a'}, 1);


% Within-cell split-half correlation
% iCell = 1;
% iMov = 1;
for iCell = 1:size(S, 1)
    tempMatRho=[];
    for iMov = 1:size(S, 2)
        if isempty(S(iCell, iMov).cellID)
            continue;
        end
        curSDF = S(iCell, iMov).matsdf{1};
        nTrial = size(curSDF, 2);
        t = 1:1:nTrial;
        oddT = 1:2:nTrial;
        setOddT = curSDF(:, ismember(t, oddT));
        setEvenT = curSDF(:, ~ismember(t, oddT));
        rho = corr(mean(setOddT, 2), mean(setEvenT, 2), 'rows', 'complete');
        tempMatRho(iMov,1) = rho;        
        
        fprintf(1, 'Computing %s : within-cell correlation for movie #%d: cell #%d/%d)..\n', ...
            nameSubjNeural, setMovie(iMov), iCell, size(S,1))
        splitHalfCorr(iMov).WithinCells(iCell).cellID = S(iCell, iMov).cellID;
        splitHalfCorr(iMov).WithinCells(iCell).movID = S(iCell, iMov).movID;
        splitHalfCorr(iMov).WithinCells(iCell).matRho = tempMatRho;
        splitHalfCorr(iMov).WithinCells(iCell).nTrial = nTrial;
    end    
end


% Across cells split-half correlation
matPair = nchoosek(1:1:size(S, 1), 2);
nPair = size(matPair, 1);
for iMov = 1:size(S, 2)
    for iPair = 1:nPair        
        if isempty(S(matPair(iPair, 1), iMov).cellID) || isempty(S(matPair(iPair, 2), iMov).cellID)
            splitHalfCorr(iMov).BetweenCells(iPair).exist = 0;
            continue;
        end
        splitHalfCorr(iMov).BetweenCells(iPair).exist = 1;
        for iCell = 1:2
            clear setOddT setEvenT
            
            idCell = matPair(iPair, iCell);
%             if isempty(S(idCell, iMov).cellID)
%                 continue;
%             end
            curSDF = S(idCell, iMov).matsdf{1};
            nTrial = size(curSDF, 2);
            t = 1:1:nTrial;
            oddT = 1:2:nTrial;
            setOddT = curSDF(:, ismember(t, oddT));
            setEvenT = curSDF(:, ~ismember(t, oddT));
            
            matSDF{iCell} = cat(2, mean(setOddT,2), mean(setEvenT,2));
        end
        fprintf(1, 'Computing %s : between-cell correlation for movie #%d: cell pair %d and %d (pair #%d / %d)..\n', ...
            nameSubjNeural, setMovie(iMov), matPair(iPair,1), matPair(iPair,2), iPair, nPair)
        rho = corr(matSDF{1}, matSDF{2}, 'rows', 'complete');
        splitHalfCorr(iMov).BetweenCells(iPair).cellID_pair = cellstr(cat(1,S(matPair(iPair,:), iMov).cellID));
        splitHalfCorr(iMov).BetweenCells(iPair).movID = S(idCell, iMov).movID;
        splitHalfCorr(iMov).BetweenCells(iPair).matRho = rho;
    end
    fprintf(1, 'Done for movie #%d! \n', setMovie(iMov))
end

% Save the results
saveFileName = fullfile(dirDataNeural, sprintf('%s_splitHalfCorr_indMov_movie%s.mat', nameSubjNeural, num2str(setMovie, '%d%d%d')));
save(saveFileName, 'splitHalfCorr', 'paramSplitHalfCorr')

% Save the figures
% Quick plot of histogram
matWC=[];
for iM = 1:length(splitHalfCorr)
tempWC = splitHalfCorr(iM).WithinCells;
matWC = cat(1, matWC, cat(1, tempWC.matRho));
end
matWC = matWC(abs(matWC)>0);

matBC = [];
for iM = 1:length(splitHalfCorr)
tempBC = splitHalfCorr(iM).BetweenCells;
indValid = cat(1, tempBC.exist)==1;
matBC = cat(1, matBC, cat(1, tempBC(indValid).matRho));
end

x = -1:0.1:1;
nWC = hist(matWC, x);
nBC = hist(matBC(:), x);

figHisto = figure;
set(figHisto, 'Color', 'w', 'PaperPositionMode', 'auto')
subplot(2, 1, 1)
bar(x, cat(2, nWC', nBC'), 'grouped')
xlim([-1.1 1.1])
legend('within', 'between', 'Location', 'Best')
title(sprintf('%s: split-half correlation', nameSubjNeural))
subplot(2, 1, 2)
bar(x, cat(2, (nWC./max(nWC))', (nBC./max(nBC))'), 'grouped')
xlim([-1.1 1.1])
legend('within', 'between', 'Location', 'Best')
title(sprintf('%s: split-half correlation: normalized', nameSubjNeural))

dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
print(figHisto, fullfile(dirFig, sprintf('splitHalfCorr_%s', nameSubjNeural)), '-depsc')



