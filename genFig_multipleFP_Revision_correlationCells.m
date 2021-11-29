% genFig_multipleFP_Revision_correlationCells.m
%
% 2021/11/1 SHP
% Compute correlations across neurons in each area using multiple temporal
% scales
% This code is generated during a revision process of Park et al., 2021 Sci Adv.

clear all;

%% Settings
flagBiowulf = 1; %0;

if flagBiowulf
    directory.dataHome = '/data/parks20/procdata/NeuroMRI/';
    dirFig = '/data/parks20/analysis/_figs';
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        directory.projects = '/Volumes/PROJECTS';
        directory.procdata = '/Volumes/PROCDATA';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/Volumes/LIBRARY';
        addpath(fullfile(directory.library, 'matlab_utils'));
    else % on virtual machine
        directory.projects = '/nifvault/NIFVAULT/projects';
        directory.procdata = '/procdata';
        directory.dataHome = fullfile(directory.procdata, 'parksh', '_macaque');
        directory.library = '/library';
        addpath(fullfile(directory.library, 'matlab_utils'));
    end
end
dirFig = fullfile(directory.projects, 'parksh/NeuroMRI/_labNote/_figs');

setNameSubjNeural = {'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan', 'Moc', 'Was', 'Dav'};

load(fullfile(directory.dataHome, sprintf('matRaster_Movie123_allCells.mat')), 'matTS_FP')
load('/procdata/parksh/_macaque/multipleFP_fsi.mat')
locFaceCell =  find(fsi.matFSI(:,1)>0.33); % find(abs(fsi.matFSI(:,1))>0.33);

matRaster = matTS_FP.matRaster(:, locFaceCell); 

catChanID = matTS_FP.catChanID(locFaceCell); % catChanID;
catSubjID = matTS_FP.catSubjID(locFaceCell); % catChanID;
setArea = matTS_FP.setArea; % catChanID;
catAreaID = matTS_FP.catAreaID(locFaceCell); % catChanID;

setSigma = [5 20 50 2400]; 

% load the TR-resolution time series 
clear matTS_FP
load(fullfile(directory.dataHome, sprintf('matSDF_Movie123_allCells.mat')), 'matTS_FP')
matFR_TR = matTS_FP.matFR_TR(:, locFaceCell); 
matNeuralRGR = matTS_FP.matNeuralRGR(:, locFaceCell); 

clear matTS_FP

%% Compute correlation coefficients using SDFs in different temporal resolution for each area 
corrSDF = struct([]); matSDF_area = struct([]);
for iArea = 1:length(setArea)
    clear cur*
    curRaster = matRaster(:, catAreaID==iArea);
    
    for iSigma = 1:length(setSigma)
        % spike density function
        curSigma  = setSigma(iSigma);   % gaussian kernel SD, in ms
        k = normpdf(-3*curSigma:3*curSigma, 0, curSigma)'; %normpdf(-500:500, 0, curSigma)'; %normpdf(-40:40, 0, curSigma)';
        
        cursdf     = conv2(curRaster, k, 'same').*1000;
        
%         matSDF{iArea, iSigma} = cursdf;
        [tempMatR] = corr(cursdf, 'type', 'spearman');
        ind = tril(true(size(tempMatR)), -1);
        curR = tempMatR(ind);
        
        corrSDF(iArea, iSigma).matR = curR;
        corrSDF(iArea, iSigma).median = median(curR);
        corrSDF(iArea, iSigma).mean = mean(curR);
        corrSDF(iArea, iSigma).ste = std(curR)./sqrt(length(curR)-1);
        corrSDF(iArea, iSigma).numCell = size(curRaster, 2);
        corrSDF(iArea, iSigma).sigma = curSigma;
        corrSDF(iArea, iSigma).catChanID = catChanID(catAreaID==iArea);
        corrSDF(iArea, iSigma).catSubjID = catSubjID(catAreaID==iArea);
        
        matSDF_area(iArea, iSigma).matSDF = cursdf;
        matSDF_area(iArea, iSigma).sigma = curSigma;
        matSDF_area(iArea, iSigma).numCell = size(curRaster, 2);
        matSDF_area(iArea, iSigma).catChanID = catChanID(catAreaID==iArea);
        matSDF_area(iArea, iSigma).catSubjID = catSubjID(catAreaID==iArea);
    end
    
%     % TR-resolution, convolved with hemodynamic function
%     
%     %     curFR_TR = matFR_TR(:, catAreaID==iArea);
%     curNeuralRGR = matNeuralRGR(:, catAreaID==iArea);
%     
%     [tempMatR] = corr(curNeuralRGR, 'type', 'spearman', 'rows', 'complete');
%     ind = tril(true(size(tempMatR)), -1);
%     curR = tempMatR(ind);
%     
%     corrSDF(iArea, length(setSigma)+1).matR = curR;
%     corrSDF(iArea, length(setSigma)+1).median = median(curR);
%     corrSDF(iArea, length(setSigma)+1).mean = mean(curR);
%     corrSDF(iArea, length(setSigma)+1).numCell = size(curRaster, 2);
%     corrSDF(iArea, length(setSigma)+1).sigma = 2400;
    
end

%% Randomly selecting pairs from different areas
% prepare indices of cells for each area
for iArea = 1:4
    indArea{iArea} = find(catAreaID==iArea);
end

nTarget = 2415; % target number of pairs (maximum number of pairs from within-area calculation)

setPairs = nchoosek(1:4, 2);
corrSDF_randAcrossArea = struct([]);
iCount = 0;
for iChoose = 1:ceil(nTarget./length(setPairs))
    for iPair = 1:length(setPairs)
        iCount = iCount+1;
        
        tL(1) = randperm(length(indArea{setPairs(iPair, 1)}), 1); % random selection of cell from area 1
        tL(2) = randperm(length(indArea{setPairs(iPair, 2)}), 1); % random selection of cell from area 2

        % raster of the randomly selected pair of neurons from this set of areas
        curRaster = matRaster(:, [indArea{setPairs(iPair, 1)}(tL(1)), indArea{setPairs(iPair, 2)}(tL(2))]);
        
        for iSigma = 1:length(setSigma)
            % spike density function
            curSigma  = setSigma(iSigma);   % gaussian kernel SD, in ms
            k = normpdf(-3*curSigma:3*curSigma, 0, curSigma)'; %normpdf(-500:500, 0, curSigma)'; %normpdf(-40:40, 0, curSigma)';
            
            cursdf     = conv2(curRaster, k, 'same').*1000;
            
            [tempR] = corr(cursdf, 'type', 'spearman');
            
            corrSDF_randAcrossArea(iSigma).matR(iCount, 1) = tempR(2,1); % add this particular correlation coefficient to the matrix
            corrSDF_randAcrossArea(iSigma).areaPair(iCount, :) = setPairs(iPair, :);
        end
    end
end

for iSigma = 1:length(setSigma)
    corrSDF_randAcrossArea(iSigma).median = median(corrSDF_randAcrossArea(iSigma).matR);
    corrSDF_randAcrossArea(iSigma).mean = mean(corrSDF_randAcrossArea(iSigma).matR);
    corrSDF_randAcrossArea(iSigma).ste = std(corrSDF_randAcrossArea(iSigma).matR)./sqrt(length(corrSDF_randAcrossArea(iSigma).matR)-1);
    corrSDF_randAcrossArea(iSigma).sigma_ms = setSigma(iSigma);
end

save('/procdata/parksh/_macaque/multipleFP_Revision_corrSDF.mat', 'corrSDF*', 'matSDF_area', '-v7.3')

    
%% Plot the distribution of correlation
setArea = {'AF', 'pAM', 'aAM', 'ML'}; 

figHistCorr = figure;
set(figHistCorr, 'Color', 'w', 'Position', [700 700 680 640]);

orderArea = [4 1 2 3];
orderSigma = [4 1 2 3];
for iArea = 1:size(corrSDF, 1)
    idArea = orderArea(iArea);
    for iSigma = 1:size(corrSDF, 2)
        idSigma = orderSigma(iSigma);
        sp(size(corrSDF, 2)*(iSigma-1)+iArea) = subplot(length(setArea), size(corrSDF, 2), size(corrSDF, 2)*(iSigma-1)+iArea);
        histogram(corrSDF(idArea, idSigma).matR, 20, 'EdgeColor', 'none');
        hold on;
        ylim = get(gca, 'YLim');
%         line([corrSDF(iArea, idSigma).median, corrSDF(iArea, idSigma).median], ylim, 'Color', 'r', 'LineWidth', 1.5);
        set(gca, 'ylim', ylim, 'YTick', ylim, 'YTickLabel', ylim);
    end
end

set(sp(:), 'XLim', [-1 1])
set(sp, 'Box', 'off', 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k')

for iArea = 1:size(corrSDF, 1)
    curSetSP = [0 4 8 12] + iArea;
    catYLim = cat(1, sp(curSetSP).YLim);
    newYlim = [0 max(catYLim(:,2))];
    set(sp(curSetSP), 'YLim', newYlim, 'YTick', newYlim, 'YTickLabel', newYlim);
    
    idArea = orderArea(iArea);
    for ii = 1:4
        idSigma = orderSigma(ii);
        plot(sp(curSetSP(ii)), corrSDF(idArea, idSigma).median, max(newYlim), 'rv', 'MarkerSize', 6, 'MarkerFaceColor', 'r');
        line(sp(curSetSP(ii)), [0 0], newYlim, 'Color', 'k', 'LineStyle', '--')
    end
end

% save
print(figHistCorr, fullfile(dirFig, 'fig_multipleFP_Revision_corrSDFs_distributions_colMLAFpAMaAM_row2400_5_20_50'), '-depsc');


%% Plot the median of correlation coefficients within each area 
% Colormap
cMap_Area = [179 226 205; 141 160 203; 252 141 98; 231 41 138]./255; % AF-pAM-aAM-ML

figMedianCorr = figure;
set(figMedianCorr, 'Color', 'w', 'Position', [700 700 345 495])

matMedianCorr = flipud(reshape(cat(1, corrSDF.median), length(setArea), length(setSigma))'); % sigma x area
matMeanCorr = flipud(reshape(cat(1, corrSDF.mean), length(setArea), length(setSigma))'); % sigma x area
matSTECorr = flipud(reshape(cat(1, corrSDF.ste), length(setArea), length(setSigma))'); % sigma x area
for iArea = 1:length(setArea)
    plot(matMeanCorr(:, iArea), '-', 'Color', cMap_Area(iArea,:), 'LineWidth', 2);
    hold on;
    line(repmat(1:length(setSigma), 2, 1),...
        [matMeanCorr(:, iArea) - matSTECorr(:, iArea), matMeanCorr(:, iArea) + matSTECorr(:, iArea)]', ...
        'Color', cMap_Area(iArea, :), 'LineWidth', 2);
    plot(matMeanCorr(:, iArea), 'o', 'Color', cMap_Area(iArea,:), 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
end
    
% p = plot(matMeanCorr);
set(gca, 'XLim', [1-0.5 size(matMeanCorr, 1)+0.5])
set(gca, 'XTick', 1:1:size(matMeanCorr, 1), 'XTickLabel', fliplr(string(setSigma)))
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 15, 'XColor', 'k', 'YColor', 'k')
% set(p, {'Color'}, mat2cell(cMap_Area, [1 1 1 1], [3]))
% set(p, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w')

% save
print(figMedianCorr, fullfile(dirFig, 'fig_multipleFP_Revision_meanCorrSDFs'), '-depsc');



%% Including results from random pairs across areas
% Plot the histogram
figHistCorr_randAcrossArea = figure;
set(figHistCorr_randAcrossArea, 'Color', 'w', 'Position', [1000 700 680 120]);

for iSigma = 1:length(corrSDF_randAcrossArea)
    sp(iSigma) = subplot(1, length(corrSDF_randAcrossArea), iSigma);
    histogram(corrSDF_randAcrossArea(iSigma).matR, 20, 'EdgeColor', 'none');
    hold on;
    ylim = get(gca, 'YLim');
    line([corrSDF_randAcrossArea(iSigma).median, corrSDF_randAcrossArea(iSigma).median], ylim, 'Color', 'r', 'LineWidth', 1.5);
    set(gca, 'ylim', ylim, 'YTick', ylim, 'YTickLabel', ylim);
end

set(sp(:), 'XLim', [-1 1])
set(sp, 'Box', 'off', 'TickDir', 'out', 'XColor', 'k', 'YColor', 'k')

% save
print(figHistCorr_randAcrossArea, fullfile(dirFig, 'fig_multipleFP_Revision_corrSDFs_randomPairsAcrossAreas_hist'), '-depsc');

%% Plot the median of correlation coefficients within each area 
% Colormap
cMap = [179 226 205; 141 160 203; 252 141 98; 231 41 138; 0 0 0]./255; % AF-pAM-aAM-ML

figMedianCorr_AcrossArea = figure;
set(figMedianCorr_AcrossArea, 'Color', 'w', 'Position', [700 700 345 495])

matMeanCorr_withAcrossArea = flipud(cat(2, reshape(cat(1, corrSDF.mean), length(setArea), length(setSigma))', cat(1, corrSDF_randAcrossArea.mean))); % sigma x area
matSTECorr_withAcrossArea = flipud(cat(2, reshape(cat(1, corrSDF.ste), length(setArea), length(setSigma))', cat(1, corrSDF_randAcrossArea.ste))); % sigma x area
for iArea = 1:length(matMeanCorr_withAcrossArea)
    plot(matMeanCorr_withAcrossArea(:, iArea), '-', 'Color', cMap(iArea,:), 'LineWidth', 2);
    hold on;
    line(repmat(1:length(setSigma), 2, 1),...
        [matMeanCorr_withAcrossArea(:, iArea) - matSTECorr_withAcrossArea(:, iArea), matMeanCorr_withAcrossArea(:, iArea) + matSTECorr_withAcrossArea(:, iArea)]', ...
        'Color', cMap(iArea, :), 'LineWidth', 2);
    plot(matMeanCorr_withAcrossArea(:, iArea), 'o', 'Color', cMap(iArea,:), 'LineWidth', 2, 'MarkerFaceColor', 'w', 'MarkerSize', 10);
end
    
% p = plot(matMeanCorr);
set(gca, 'XLim', [1-0.5 size(matMeanCorr_withAcrossArea, 1)+0.5])
set(gca, 'XTick', 1:1:size(matMeanCorr_withAcrossArea, 1), 'XTickLabel', string(setSigma))
set(gca, 'Box', 'off', 'TickDir', 'out', 'FontSize', 15, 'XColor', 'k', 'YColor', 'k')
% set(p, {'Color'}, mat2cell(cMap_Area, [1 1 1 1], [3]))
% set(p, 'LineWidth', 2, 'Marker', 'o', 'MarkerSize', 10, 'MarkerFaceColor', 'w')

% save
print(figMedianCorr_AcrossArea, fullfile(dirFig, 'fig_multipleFP_Revision_corrSDFs_randomPairsAcrossAreas_mean'), '-depsc');

%% number of cells
numCells = cat(1, corrSDF(:, 1).numCell);

%% Correlation across all the cells
curMatSDF = cat(2, matSDF_area(:, 4).matSDF); % 5ms window
[matR_all] = corr(curMatSDF, 'type', 'Spearman');

curMatSDF_mixed = []; % mix across monkeys
for iArea = 1:size(matSDF_area, 1)
    curMatSDF_mixed = cat(2, curMatSDF_mixed, matSDF_area(iArea, 4).matSDF(:, randperm(numCells(iArea))));
end
[matR_all_mixed] = corr(curMatSDF_mixed, 'type', 'Spearman');

figure
imagesc(matR_all_mixed)
colormap(jet)
set(gca, 'CLim', [-1 1].*0.4);
set(gca, 'XTick', cumsum(numCells), 'YTick', cumsum(numCells))

%% Correlation between single cells and population activity (averaged signal from each patch)
load(fullfile(directory.dataHome, '/multipleFP_Revision_corrSDF.mat'), 'matSDF_area')
% Let's start with ML and aAM
numCells = cat(1, matSDF_area(:, 1).numCell);

orderArea = [4 1 2 3];
for iArea = 1:4 %size(matSDF_area, 1)
    matSDF_area_5ms{iArea} = matSDF_area(orderArea(iArea), 1).matSDF;
    matSDF_area_2400ms{iArea} = matSDF_area(orderArea(iArea), 4).matSDF;
    avgSDF_area(:, iArea) = mean(matSDF_area_5ms{iArea}, 2);     
    avgSDF_area_2400ms(:, iArea) = mean(matSDF_area_2400ms{iArea}, 2);  
end

for iArea = 1:4
    R_MUA_5ms{iArea} = corr(matSDF_area_5ms{iArea}, avgSDF_area(:, iArea), 'type', 'Spearman');
    R_MUA_2400ms{iArea} = corr(matSDF_area_2400ms{iArea}, avgSDF_area_2400ms(:, iArea), 'type', 'Spearman');
end

setMedian(1,:) = cellfun(@median, R_MUA_5ms);
setMedian(2,:) = cellfun(@median, R_MUA_2400ms);

edges = -1:0.1:1;
figure;
set(gcf, 'Position', [500 720 1190 260])
for iArea = 1:4
    sss(iArea) = subplot(1, 4, iArea);
    histogram(R_MUA_5ms{iArea}, edges, 'faceColor', 'b', 'edgecolor', 'none'); hold on;
    histogram(R_MUA_2400ms{iArea}, edges, 'faceColor', 'm', 'edgecolor', 'none');    
    hold on
    line([0 0], get(gca, 'YLim'), 'Color', 'k', 'LineStyle', '--')
%     P1(iArea) = plot(setMedian(1, iArea), max(get(gca, 'YLim')), 'bv', 'MarkerFaceColor', 'b');
%     P2(iArea) = plot(setMedian(2, iArea), max(get(gca, 'YLim')), 'mv', 'MarkerFaceColor', 'm');    
end
set(sss, 'XColor', 'k', 'YColor', 'k', 'TickDir', 'out', 'Box', 'off')
print(gcf, fullfile(dirFig, 'multipleFP_Revision_corrCellsMUA_MLAFpAMaAM_distMedian'), '-depsc')

figure;
set(gcf, 'Position', [500 720 1190 260])
for iArea = 1:4
    ssss(iArea) = subplot(1, 4, iArea);
    plot(R_MUA_2400ms{iArea}, R_MUA_5ms{iArea}, 'bo', 'MarkerFaceColor', 'w', 'LineWidth', 1);
    set(gca, 'XLim', [-0.5 1], 'YLim', [-0.5 1]);
    line([0 0; -0.5 1; -0.5 1]', [-0.5 1; 0 0; -0.5 1]', 'Color', 'k', 'LIneStyle', '--') 
%     histogram(R_MUA_5ms{iArea}, edges, 'faceColor', 'b'); hold on;
%     histogram(R_MUA_2400ms{iArea}, edges, 'faceColor', 'm');       
end
set(ssss, 'XColor', 'k', 'YColor', 'k', 'TickDir', 'out', 'Box', 'off')
print(gcf, fullfile(dirFig, 'multipleFP_Revision_corrCellsMUA_MLAFpAMaAM_2400vs5'), '-depsc')

%% Cross-correlation 
load(fullfile(directory.dataHome, '/multipleFP_Revision_corrSDF.mat'), 'matSDF_area')
% Let's start with ML and aAM
numCells = cat(1, matSDF_area(:, 1).numCell);

% prepare highpass-filtered ts from each area
Fs = 1000; %1kHz sampling rate
fpass = 5; %10; %filtering criterion in Hz.

orderArea = [4 1 2 3];
for iArea = 1:4 %size(matSDF_area, 1)
    matSDF_area_5ms{iArea} = matSDF_area(orderArea(iArea), 1).matSDF;
    matSDF_area_filtered{iArea} = highpass(matSDF_area(orderArea(iArea), 1).matSDF, fpass, Fs);
    avgSDF_area(:, iArea) = mean(matSDF_area_5ms{iArea}, 2);     
end
avgSDF_area_filtered = highpass(avgSDF_area, fpass, Fs);

clear matSDF_area

% 1. Across single cell pairs from different areas
setAreaName = {'ML', 'AF', 'pAM', 'aAM'};
setAreaPair = nchoosek(1:4, 2);
maxlag = 200;
pairXCF = struct([]);
for iPair = 2:length(setAreaPair)
    indArea1 = setAreaPair(iPair, 1);
    indArea2 = setAreaPair(iPair, 2);
    
    matMaxXCF_r = NaN(size(matSDF_area_filtered{indArea2}, 2), size(matSDF_area_filtered{indArea1}, 2));
    matMaxXCF_lag = NaN(size(matSDF_area_filtered{indArea2}, 2), size(matSDF_area_filtered{indArea1}, 2));
    
    matMaxXCF_r_f = NaN(size(matSDF_area_filtered{indArea2}, 2), size(matSDF_area_filtered{indArea1}, 2));
    matMaxXCF_lag_f = NaN(size(matSDF_area_filtered{indArea2}, 2), size(matSDF_area_filtered{indArea1}, 2));
    matXR = NaN(maxlag*2+1, size(matSDF_area_filtered{indArea2}, 2), size(matSDF_area_filtered{indArea1}, 2));
    matXR_f = NaN(maxlag*2+1, size(matSDF_area_filtered{indArea2}, 2), size(matSDF_area_filtered{indArea1}, 2));

    for iC1 = 1:size(matSDF_area_filtered{indArea1}, 2)        
        curmatXR = []; curmatXR_f=[];
        for iC2 = 1:size(matSDF_area_filtered{indArea2}, 2)
            fprintf(1, 'Area %s x %s: Cell %d/%d & Cell %d/%d \n', setAreaName{indArea1}, setAreaName{indArea2}, ...
            iC1, size(matSDF_area_filtered{indArea1}, 2), iC2, size(matSDF_area_filtered{indArea2}, 2))
            
            [xcf,lags,bounds] = crosscorr(matSDF_area_5ms{indArea1}(:, iC1), matSDF_area_5ms{indArea2}(:, iC2), ...
                'NumLags', maxlag);
            [xcf_f,lags,bounds] = crosscorr(matSDF_area_filtered{indArea1}(:, iC1), matSDF_area_filtered{indArea2}(:, iC2), ...
                'NumLags', maxlag);
%             [R, lags] = xcorr(matSDF_area_filtered{indArea1}(:, iC1), matSDF_area_filtered{indArea2}(:, iC2), ...
%                 maxlag, 'normalized');
            
            curmatXR = cat(2, curmatXR, xcf);
            [maxXCF_r, maxXCF_lag] = max(xcf);
            matMaxXCF_r(iC2, iC1) = maxXCF_r; 
            matMaxXCF_lag(iC2, iC1) = lags(maxXCF_lag);
            
            curmatXR_f = cat(2, curmatXR_f, xcf_f);
            [maxXCF_r_f, maxXCF_lag_f] = max(xcf_f);
            matMaxXCF_r_f(iC2, iC1) = maxXCF_r_f;
            matMaxXCF_lag_f(iC2, iC1) = lags(maxXCF_lag_f);
        end
        matXR(:, :, iC1) = curmatXR;
        matXR_f(:, :, iC1) = curmatXR_f;
    end
    pairXCF(iPair).matXCF = matXR;
    pairXCF(iPair).matXCF_f = matXR_f;
    pairXCF(iPair).matMaxXCF_r = matMaxXCF_r; % area 2 (row) x area 1 (col)
    pairXCF(iPair).matMaxXCF_lag = matMaxXCF_lag;
    pairXCF(iPair).matMaxXCF_r_f = matMaxXCF_r_f;
    pairXCF(iPair).matMaxXCF_lag_f = matMaxXCF_lag_f;
end

% Averaged across cells within each area
setAreaName = {'ML', 'AF', 'pAM', 'aAM'};
setAreaPair = nchoosek(1:4, 2);
maxlag = 200;
for iPair = 1:length(setAreaPair)
    [xcf,lags,bounds] = crosscorr(avgSDF_area(:, setAreaPair(iPair, 1)), avgSDF_area(:, setAreaPair(iPair, 2)), ...
        'NumLags', maxlag);
    [xcf_f,lags,bounds] = crosscorr(avgSDF_area_filtered(:, setAreaPair(iPair, 1)), avgSDF_area_filtered(:, setAreaPair(iPair, 2)), ...
        'NumLags', maxlag);
    
    [maxXCF_r, maxXCF_lag] = max(xcf);    
    [maxXCF_r_f, maxXCF_lag_f] = max(xcf_f);
    
    pairXCF(iPair).xcfAvgArea = xcf;
    pairXCF(iPair).xcfAvgArea_maxXCF_r = maxXCF_r; 
    pairXCF(iPair).xcfAvgArea_maxXCF_lag = lags(maxXCF_lag);
    pairXCF(iPair).xcfAvgArea_f = xcf_f;
    pairXCF(iPair).xcfAvgArea_maxXCF_r_f = maxXCF_r_f; 
    pairXCF(iPair).xcfAvgArea_maxXCF_lag_f = lags(maxXCF_lag_f);
    pairXCF(iPair).setPair = setAreaPair(iPair, :);
    pairXCF(iPair).setAreaName = setAreaName(setAreaPair(iPair,:));
end

save('/procdata/parksh/_macaque/multipleFP_Revision_crosscorrSDF.mat', 'pairXCF')


% plot the results
figure;
for iPair = 1:length(setAreaPair)
    sp(iPair) = subplot(2, 3, iPair);
    plot(pairXCF(iPair).matMaxXCF_lag_f(:), pairXCF(iPair).matMaxXCF_r_f(:), 'bo')
    hold on
    plot(pairXCF(iPair).matMaxXCF_lag(:), pairXCF(iPair).matMaxXCF_r(:), 'ro')
    title(pairXCF(iPair).setAreaName)
end

figure;
for iPair = 1:length(setAreaPair)
    sp(iPair) = subplot(2, 3, iPair);
    histogram(pairXCF(iP).matMaxXCF_lag(:), 20); hold on;
    histogram(pairXCF(iP).matMaxXCF_lag_f(:), 20); 
    title(pairXCF(iPair).setAreaName)
end

% plot the results
figure;
for iPair = 1:length(setAreaPair)
    sp(iPair) = subplot(2, 3, iPair);
    plot(lags, pairXCF(iPair).xcfAvgArea);
    title(pairXCF(iPair).setAreaName)
end

% plot the results
figure;
for iPair = 1:length(setAreaPair)
    sp(iPair) = subplot(2, 3, iPair);
    plot(lags, pairXCF(iPair).xcfAvgArea_f, 'r-');
    title(pairXCF(iPair).setAreaName)
end


            
% figXCorr = figure;
% set(figXCorr, 'Color', 'w', 'Position', [1000 700 1000 800]);
% 
% for area1 = 1:4
%     for area2 = 1:4
%         ssp(4*(area1-1)+area2) = subplot(4, 4, 4*(area1-1)+area2); 
%         plot(lags, R(:,4*(area1-1)+area2));
%     end
% end


%% Maybe cells in each cluster? have different relationship?
%e.g. from same area but fall into different groups based on whole-brain
%network. could they have different relationship?




        
% avgSDF_area(:, 1) = mean(matSDF_area(4, 1).matSDF, 2); %ML
% avgSDF_area(:, 2) = mean(matSDF_area(1, 1).matSDF, 2); %AF
% avgSDF_area(:, 3) = mean(matSDF_area(2, 1).matSDF, 2); %pAM
% avgSDF_area(:, 4) = mean(matSDF_area(3, 1).matSDF, 2); %aAM
% 
% % highpass filtering
% % curMatSDF = cat(2, matSDF_area(:, 1).matSDF);
% Fs = 1000; %1kHz sampling rate
% fpass = 5; %10; %filtering criterion in Hz.
% avgSDF_area_filtered = highpass(avgSDF_area, fpass, Fs);
% 
% % Use crosscorr function for each pair
% setAreaPair = nchoosek(1:4, 2);
% maxlag = 300;  %maximum lag in ms
% setXCF = NaN(maxlag*2+1, size(setAreaPair, 1));
% for iPair = 1:length(setAreaPair)
%     indArea1 = setAreaPair(iPair, 1);
%     indArea2 = setAreaPair(iPair, 2);
%     
%     [xcf,lags,bounds] = crosscorr(avgSDF_area_filtered(:, indArea1), avgSDF_area_filtered(:, indArea2), 'NumLags', maxlag);
%     setXCF(:, iPair) = xcf;
% end


% % Use "xcorr" function
% maxlag = 100; %300;
% [R, lags] = xcorr(avgSDF_area_filtered, maxlag, 'normalized');
% 
% figXCorr = figure;
% set(figXCorr, 'Color', 'w', 'Position', [1000 700 1000 800]);
% 
% for area1 = 1:4
%     for area2 = 1:4
%         ssp(4*(area1-1)+area2) = subplot(4, 4, 4*(area1-1)+area2); 
%         plot(lags, R(:,4*(area1-1)+area2));
%     end
% end
        
% Okay, what about when face appears? maybe looking at only 1-s after big face
% appears?
s = []; e = []; 
matXR_face = [];
for iMovie = 1:3
load(sprintf('/procdata/parksh/MovieRegressors/annotationMovie%d.mat', iMovie)) % DM's scene annotation

faceA = [];
for iS = 1:length(epoch)
    faceA(iS,1) = epoch(iS).notes.face.A; % face area
end
tLoc = find(faceA > 30); % arbitrary criterion of 10 to select scenes with larger face
s = cat(1, s, (300*30*(iMovie-1)+sta(tLoc)').*33 - 30); % start slightlly early by rounding off the 0.33 ms from 33.33ms/frame
e = cat(1, e, s + 2000); % end of 1-s window in ms

% L = length(s);
% F = cumsum(e-s+1);
% idx = ones(1, F(end));
% idx(1) = s(1);
% idx(1+F(1:L-1)) = s(2:L)-e(1:L-1);
% idx = cumsum(idx);

for iSS = 1:length(tLoc)
    maxlag = 300;
%     [xcf,lags,bounds] = crosscorr(avgSDF_area(:, indArea1), avgSDF_area(:, indArea2), 'NumLags', maxlag);
    [R, lags] = xcorr(avgSDF_area_filtered(s(iSS):e(iSS), :), maxlag, 'normalized');
    matXR_face = cat(3, matXR_face, R);
end

end

meanXR = mean(matXR_face, 3); %matXR_face(:, :, 2); %mean(matXR_face, 3);

figXCorr = figure;
set(figXCorr, 'Color', 'w', 'Position', [1000 700 1000 800]);
size
for area1 = 1:4
    for area2 = 1:4
        ssp(4*(area1-1)+area2) = subplot(4, 4, 4*(area1-1)+area2); 
        plot(lags, meanXR(:,4*(area1-1)+area2));
    end
end









