% genFig_SFN2018.m
%
% 2018/10/15 SHP
% 1) Trying out various pilot analyses and 2) making figures for the SFN 2018

%% Plot of cell correlation
% Save the results
setNameSubj = {'Dav', 'Tor', 'Spi', 'Dan', 'Mat'}; % ML AF AM

for iS = 1:length(setNameSubj)
    nameSubjNeural = setNameSubj{iS};
    dirDataNeural = fullfile('/procdata/parksh', nameSubjNeural);
    setMovie = [1 2 3 4 5 6];
    if sum(strcmpi(nameSubjNeural, {'tor', 'toroid'}))
        setMovie = [1 2 3];
    end
    load(fullfile(dirDataNeural, sprintf('%s_splitHalfCorr_indMov_movie%s.mat', nameSubjNeural, num2str(setMovie, '%d%d%d'))))
    
    matWC=[];
    for iM = 1:length(splitHalfCorr)
        clear tempWC
        tempWC = splitHalfCorr(iM).WithinCells;
        matWC = cat(1, matWC, cat(1, tempWC.matRho));
    end
    matWC = matWC(abs(matWC)>0);
    
    matBC = [];
    for iM = 1:length(splitHalfCorr)
        clear tempBC
        tempBC = splitHalfCorr(iM).BetweenCells;
        indValid = cat(1, tempBC.exist)==1;
        matBC = cat(1, matBC, cat(1, tempBC(indValid).matRho));
    end
    
    x = -1:0.1:1;
    nWC = hist(matWC, x);
    nBC = hist(matBC(:), x);
    
    indSubj(iS).name = nameSubjNeural;
    indSubj(iS).histX = x;
    indSubj(iS).nWC = nWC;
    indSubj(iS).nBC = nBC;
    indSubj(iS).matWC = matWC;
    indSubj(iS).matBC = matBC;
end

% merge for each area
catNWC = cat(1, indSubj(1).nWC, sum(cat(1, indSubj(2:3).nWC)), sum(cat(1, indSubj(4:5).nWC)))'; % 1 (Dav) for ML, 2:3 (Tor, Spi) for AF, 4:5 (Dan, Mat) for AM
catNBC = cat(1, indSubj(1).nBC, sum(cat(1, indSubj(2:3).nBC)), sum(cat(1, indSubj(4:5).nBC)))'; % 1 (Dav) for ML, 2:3 (Tor, Spi) for AF, 4:5 (Dan, Mat) for AM
pos = [100 100 780 160];

figHisto_within = figure;
set(figHisto_within, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', pos)
x = -1:0.1:1;
b = bar(x, catNWC./repmat(max(catNWC), size(catNWC, 1), 1), 1.5);
hold on
 line([0 0], [0 1], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
xlim([-1.1 1.1]); ylim([0 1.1])
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XTick', -1:0.2:1, 'XTickLabel', [])

dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
print(figHisto_within, fullfile(dirFig, 'MLAFAM_splitHalfCorr_withinCell_normFreq'), '-depsc')


figHisto_between= figure;
set(figHisto_between, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', pos)
x = -1:0.1:1;
bb = bar(x, catNBC./repmat(max(catNBC), size(catNBC, 1), 1), 1.5);
xlim([-1.1 1.1]); ylim([0 1.1])
hold on
 line([0 0], [0 1], 'LineStyle', ':', 'LineWidth', 2, 'Color', 'k')
set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XTick', -1:0.2:1, 'XTickLabel', [])

dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
print(figHisto_between, fullfile(dirFig, 'MLAFAM_splitHalfCorr_betwenCell_normFreq'), '-depsc')


%% Dexter's pulvinar cells SDF from each trials
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
nameSubjNeural = 'Dex';% 'Mat';
% Load data files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
load(fullfile(dirDataNeural, filenameNeural), 'paramSDF') % to get the cell ID
setCellID = paramSDF.setCellIDs;
setMovie = 1;
sigma = 1000;
S = createCellRegressor_indMov(dirDataNeural, setCellID, setMovie, sigma); %{'065a', '082a', '089a'}, 1);
tempC = bone(6); % jet(size(S(1).matsdf{1}, 2))*0.75; % temporary ColorOrder for line plot

S_1 = S(13); % '20151118020a'
S_2 = S(21); % '20160414007a' % should omit the last trial

figSDF = figure;
set(figSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 600])
% case of two example neurons
SP(1) = subplot('position', [0.1 0.55 0.8 0.35]);
SP(2) = subplot('position', [0.1 0.15 0.8 0.35]);
set(SP, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(SP(1), S_1.matsdf{1}, 'LineWidth', 2)
plot(SP(2), S_2.matsdf{1}(:, 1:4), 'LineWidth', 2) %omit the last trial
set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 20,  'XTickLabel', [])

setCellIDs = {'20151118020a', '20160414007a'};
print(figSDF, fullfile(dirFig, sprintf('exampleSDF_%s_cell%s_trialbytrial', nameSubjNeural, cell2mat(setCellIDs))), '-depsc')
    
%% Davida
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
nameSubjNeural = 'Dan'; %'Dav';% 'Mat';
% Load data files
dirDataHome = '/procdata/parksh/';
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
% filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% load(fullfile(dirDataNeural, filenameNeural), 'paramSDF') % to get the cell ID
setCellID = {'18', '28'};% {'07', '17'}; %paramSDF.setCellIDs;
setMovie = 1; %4; % 1;
sigma = 1000;
S = createCellRegressor_indMov(dirDataNeural, setCellID, setMovie, sigma); %{'065a', '082a', '089a'}, 1);
tempC = bone(15); % jet(size(S(1).matsdf{1}, 2))*0.75; % temporary ColorOrder for line plot

S_1 = S(1); % '07'
S_2 = S(2); % '17' % should omit the last trial

figSDF = figure;
set(figSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 600])
% case of two example neurons
SP(1) = subplot('position', [0.1 0.55 0.8 0.35]);
SP(2) = subplot('position', [0.1 0.15 0.8 0.35]);
set(SP, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(SP(1), S_1.matsdf{1}, 'LineWidth', 2)
plot(SP(2), S_2.matsdf{1}, 'LineWidth', 2) %omit the last trial
set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 20,  'XTickLabel', [])

print(figSDF, fullfile(dirFig, sprintf('exampleSDF_%s_cell%s_trialbytrial', nameSubjNeural, cell2mat(setCellID))), '-depsc')



%
figSDF = figure;
set(figSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 300]); %
% case of five example neurons
startY = 0.09;
height = 0.24;
gap = 0.07;
for iCell = 1:length(setCellID)
    figure(figSDF); clf;
    
    set(gca, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
    plot(S(iCell, setMovie(iMovie)).matsdf{1}, 'LineWidth', 2)
    title(sprintf('Movie #%d, Cell #%d: %s', setMovie(iMovie), iCell, char(setCellID{iCell})))
    set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
    set(gca, 'XTickLabel', 0:50:300)
    xlabel('Time (s)') 
    ylabel(gca, 'Spikes / s')
    
        input('')
end



%     
%     set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
%     set(SP(2:3), 'XTickLabel', [])
%     set(SP(1), 'XTickLabel', 0:50:300)
%     xlabel(SP(1), 'Time (sec)')
%     ylabel(SP(2), 'Spike / sec')
    


%     % save as eps
%     print(figSDF, sprintf(fullfile(dirFig, '%s_Cell%sMovie%s_SDF'), nameSubjNeural, char(setCellID{iCell}), num2str(setMovie, '%d%d%d')), '-depsc')
    
% end

figSDF = figure;
set(figSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 600])
% case of two example neurons
SP(1) = subplot('position', [0.1 0.55 0.8 0.35]);
SP(2) = subplot('position', [0.1 0.15 0.8 0.35]);
set(SP, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(SP(1), S_1.matsdf{1}, 'LineWidth', 2)
plot(SP(2), S_2.matsdf{1}, 'LineWidth', 2)
set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)


%% PCA on fMRI data
% load the masks (movie-driven & brain-only mask)
load('/procdata/parksh/Art/Art_MaskArrays.mat', 'movieDrivenAmp', 'brainMask_BlockAna3D');
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D

load('/procdata/parksh/Art/Art_movieTS_fMRI_indMov.mat')

indMovieBOLD = [1 2 3];

for iM = 1:3;
    fmritc=[];
    curvoltc = voltcIndMov{iM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc = pcvoltc(:,:,:,8:125); %pcvoltc;
    
    [nx, ny, nz, nt] = size(fmritc);
    nVox = nx*ny*nz;
    
    matTS = reshape(fmritc, nVox, nt);
    matTS_moviemask = matTS(moviemask_vec,:); %matR_SU(moviemask_vec,:); % 15495 voxels
    
    [coeff, score, latent, tsquared, explained] = pca(zscore(matTS_moviemask), 'Economy', 'off', 'Centered', 'on');
    [residuals] = pcares(zscore(matTS_moviemask), 1);
    
    resultsPCA(iM).coeff = coeff(:,1:10);
    resultsPCA(iM).score = score(:,1:10);
    resultsPCA(iM).explained = explained;
    paramPCA.option = {'Economy', 'off', 'Centered', 'on'};
    paramPCA.flag_zscore = 1;
    
    resultsPCAres(iM).residuals = residuals;
    paramPCAres.ndim = 1;
    paramPCAres.flag_zscore = 1;
end

% concatenate the first principal component across movies
catCoeff = cat(1, resultsPCA.coeff);
pc1_fMRI = [];
for iM = 1:3
    pc1_fMRI = cat(1, pc1_fMRI, NaN(7,1), catCoeff(118*(iM-1)+1:118*iM, 1));
end
% pc1_fMRI = catCoeff(:,1);

save('/procdata/parksh/Art/Art_movieTS_fMRI_Movie123_PCA.mat', 'paramPCA', 'resultsPCA', 'paramPCAres', 'resultsPCAres', 'pc1_fMRI')

clear fmritc pcvoltc avgvoltc curvoltc voltcIndMov

% some evaluation on PCs
catVar = cat(2, resultsPCA.explained);
cumsum(catVar(1:10, :))



%% Correlation between neural time series and 1st PC of fMRI data
% get the neuronal time series
setMovieset = [1 2 3; 4 5 6];
setFileName = {'matSDF_DavSpiMatDan_Movie123.mat', 'matSDF_AvaDavSpiMatDan_Movie456.mat'};
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

iMovieSet = 1; %2; %1;
setMovie = setMovieset(iMovieSet,:);
tempS = num2str(setMovie);
MovieStr = tempS(~isspace(tempS));
load(fullfile('/procdata/parksh/Dan', setFileName{iMovieSet}))

% for movie 123
% for Spice: excChanInd = [6 11 30 35]; excChanID = {'09, '14', '42', '47'}
% % Cells that should be excluded indvalid = setdiff(1:40, excChanInd);
% for Davida: indvalid = [1 4 5 8 9 10]; validchanID = paramCorr.validChanID(indvalid,:); %{'02' '08' '09' '17' '18' '20'}
% for Dango:  indvalid = [3 4 5 6 7 8 9];  validchanID =paramCorr.validChanID(indvalid,:) %{'10' '13' '15' '17' '18' '23' '26'}

% for movie 456
% for Spice: excChanInd = [4 7 8 9 14 15 18 23 32 34 38 41 42 43]; indvalid = setdiff(1:49, excChanInd);
% for Davida: indvalid = [2 3 5 9 10 12 14]; validchanID = paramCorr.validChanID(indvalid,:) %  {'02', '03', '08', '17', '18', 20', '22'};
% for Dango: indvalid = [1 4 5 6 7 8]; validchanID = paramCorr.validChanID(indvalid,:) % {'04', '09', '10', '17', '23', '26'};

setIndArea = [1 1 2 3 3]; %[0 1 2 3 3]; %for matSDF_DavSpiMatDan_Movie123.mat % [1 1 2 3 3]; %for matSDF_AvaDavSpiMatDan_Movie456.mat % 1 for ML, 2 for AF, 3 for AM
catIndArea=[];
catIndSubj = [];
catCellID=[];

matFR_SU_norm =[];
matFR_SU=[];
matNeuralRGR =[];
for iSubj = 2:5 %1:5
    
    switch lower(matSDF(iSubj).nameSubj)
        case 'spi'
            if iMovieSet == 1
                excChanIndex = [6 11 30 35];
            elseif iMovieSet == 2
                excChanIndex = [4 7 8 9 14 15 18 23 32 34 38 41 42 43];
            end
            indValid = setdiff(1:length(matSDF(iSubj).setCellIDs), excChanIndex)';
        case 'dav'
            if iMovieSet == 1
                indValid = [1 4 5 8 9 10];
            elseif iMovieSet ==2
                indValid = [2 3 5 9 10 12 14];
            end
        case 'dan'
            if iMovieSet == 1
                indValid = [3 4 5 6 7 8 9];
            elseif iMovieSet == 2
                indValid = [1 4 5 6 7 8];
            end
        otherwise
            indValid = 1:length(matSDF(iSubj).setCellIDs);
    end
    
    matFR_SU_norm = cat(2, matFR_SU_norm, matSDF(iSubj).matFR_SU_norm(:, indValid));
    matNeuralRGR = cat(2, matNeuralRGR, matSDF(iSubj).matNeuralRGR(:, indValid));
    matFR_SU = cat(2, matFR_SU, matSDF(iSubj).matFR_SU(:, indValid));
    
    catIndArea = cat(1, catIndArea, ones(length(indValid), 1).*setIndArea(iSubj));
    catIndSubj = cat(1, catIndSubj, ones(length(indValid), 1).*iSubj);
    catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', matSDF(iSubj).area, matSDF(iSubj).nameSubj), length(indValid), 1)), matSDF(iSubj).setCellIDs(indValid)));
        
    setSubjID = {'Ava', 'Dav', 'Spi', 'Mat', 'Dan'};
end

flagMerge = 1;
if iMovieSet == 1 && flagMerge
    load('/procdata/parksh/Spi/2016Nov_movie/matSDF_TorRhoSigSpi_Movie123_new.mat') %% Merge earlier AF datamatSDF_TorRhoSigSpiMovie123.mat') %% Merge earlier AF data
    % need to check later about Spice cells here (whether I need to exclude
    % some cells here or not)
    
    % Some info of the earlier dataset
    setSubj_2016 = {'Tor', 'Rho', 'Sig', 'Spi'};
    setIndSubj_2016 = [6 7 8 9];
    for iS = 1:length(matSDF)
        nCell_2016(iS) = length(matSDF(iS).setCellIDs);
        
        catIndArea = cat(1, catIndArea, ones(length(matSDF(iS).setCellIDs), 1)*2); % 2 for AF (1 for ML, 3 for AM)
        catIndSubj = cat(1, catIndSubj, ones(length(matSDF(iS).setCellIDs), 1).*setIndSubj_2016(iS));
        catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', matSDF(iS).area, matSDF(iS).nameSubj),  length(matSDF(iS).setCellIDs), 1)), matSDF(iS).setCellIDs));
%         catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', 'AF', setSubj_2016{iS}), length(matSDF(iS).setCellIDs), 1)), matSDF(iS).setCellIDs));
    end
    
    matFR_SU_norm = cat(2, matFR_SU_norm, cat(2, matSDF.matFR_SU_norm));    
    matFR_SU = cat(2, matFR_SU, cat(2, matSDF.matFR_SU));    
    matNeuralRGR = cat(2, matNeuralRGR, cat(2, matSDF.matNeuralRGR));
    setSubjID = cat(2, setSubjID, setSubj_2016);
end

% Just to check: interim PCA on properly standardized neuronal data
[coeff, score, latent, tsquared, explained] = pca(zscore(matFR_SU_norm'), 'Economy', 'off', 'Centered', 'on');

%% Correlation between neural regressors and the 1st PC of the fMRI data
load('/procdata/parksh/Art/Art_movieTS_fMRI_Movie123_PCA.mat', 'pc1_fMRI')

[matR] = corr(pc1_fMRI, matNeuralRGR, 'rows', 'complete', 'type', 'Spearman');
% [Rvals, Pvals] = corr(reshape(fmritc, nVox, nt)', neuralrgrs',...
%         'rows','complete', 'type', 'Spearman');

