% genFigs_multiplePatchDataMining.m
%
% 2018/09/05 SHP
% Cells from multiple patches data mining

%% Settings
flagBiowulf = 0;

if flagBiowulf
    dirDataHome = '/data/parks20/procdata/NeuroMRI/';
    addpath('/data/parks20/analysis/NeuroMRI/'); % to use doConv.m function
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/PROJECTS';
        dirProcdata = '/Volumes/PROCDATA';
        dirDataHome = fullfile(dirProcdata, 'parksh');
        dirLibrary = '/Volumes/LIBRARY';
        addpath(fullfile(dirLibrary, 'matlab_utils'));
    else % on virtual machine
        dirProjects = '/projects';
        dirProcdata = '/procdata';
        dirDataHome = fullfile(dirProcdata, 'parksh');
        dirLibrary = '/library';
        addpath(fullfile(dirLibrary, 'matlab_utils'));
    end
end

%% Prepare time-series matrices in different resolution
setNameSubjNeural = {'Ava', 'Dav', 'Spi', 'Mat', 'Dan'}; % {'Tor', 'Rho', 'Sig', 'Spi'};
setNameArea = {'ML', 'ML', 'AF', 'AM', 'AM'};
setMovie = [1 2 3]; %[4 5 6]; %[1 2 3];

matFR_dT = struct([]);
for iSubj = 1:length(setNameSubjNeural)
    
    % Load the data
    nameSubjNeural = setNameSubjNeural{iSubj};
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
    load(fullfile(dirDataNeural, filenameNeural))
    fprintf(1, '\nLoading single unit data of %s: %s ....\n ', nameSubjNeural, filenameNeural)
        
    [indDataMat, CellID, movieID] = genDataMatrix_SU(nameSubjNeural, 0); % data matrix
    validC = find(indDataMat*ismember(movieID, setMovie)>0); % valid channel with movie [1 2 3]
    
    if isempty(validC)
        continue;
    end
    
    matFR_dT(iSubj).nameSubj = nameSubjNeural;
    matFR_dT(iSubj).setCellIDs = paramSDF.setCellIDs(validC);
    
    indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
    
    % 3. Fine temporal resolution
    dirDataNeural_individualCellFile = dirDataNeural;
    if sum(strcmpi(nameSubjNeural, {'spice', 'spi'}))
        dirDataNeural_individualCellFile = fullfile(dirDataNeural, '2018Jan_movie');
    end
    FR_dTfine = createCellRegressor_indMov_discreteTime(dirDataNeural_individualCellFile, paramSDF.setCellIDs(validC), ... %cellstr(paramCorr.validChanID),...
        setMovie, 0.1); % in 10Hz (number of spikes for every 100ms)
    % concatenate across movies
    matSU_mnFR=[];
    for iUnit = 1:size(FR_dTfine,1)
        tempFR = cat(1, FR_dTfine(iUnit, :).mnFR);
        matSU_mnFR(:,iUnit) = tempFR;
        matFR_dT(iSubj).matSU_matFR(:,iUnit)  = cat(1, FR_dTfine(iUnit,:).matFR);
    end
    matFR_dT(iSubj).matSU_mnFR = matSU_mnFR;
%     matFR_dT(iSubj).matSU_matFR = matSU_matFR;
end

%% Across-trial variance in each individual neurons: using Coefficient of variation
% iSubj=4;
% % tempVar = cellfun(@var, matFR_dT(iSubj).matSU_matFR, num2cell(zeros(size(matFR_dT(iSubj).matSU_matFR))), ...
% %     num2cell(ones(size(matFR_dT(iSubj).matSU_matFR)).*2), 'Uniformoutput', false);
% % tempMn = cellfun(@mean, tempVar);
% mnFR_nonZero = matFR_dT(iSubj).matSU_mnFR;
% mnFR_nonZero(matFR_dT(iSubj).matSU_mnFR==0) = realmin;
% 
% tempStd = cellfun(@std, matFR_dT(iSubj).matSU_matFR, num2cell(zeros(size(matFR_dT(iSubj).matSU_matFR))), ...
%     num2cell(ones(size(matFR_dT(iSubj).matSU_matFR)).*2), 'UniformOutput', false);
% tempCV = cellfun(@rdivide, tempStd, mat2cell(mnFR_nonZero, [3000 3000 3000]', ones(1, size(tempStd, 2))), 'UniformOutput', false);
% [y, ind]=sort(cat(1, tempStd{:,2})); ind(end)
% aaa = cell2mat(tempCV);
% figure
% imagesc(aaa')
% title(sprintf('%s', matFR_dT(iSubj).nameSubj))
% 
% figure
% hist(aaa(:))
% title(sprintf('%s', matFR_dT(iSubj).nameSubj))
% 
% [y, order]=sort(nanmean(aaa), 'ascend');
% cat(1, matFR_dT(iSubj).setCellIDs{order(1:5)})
% cat(1, matFR_dT(iSubj).setCellIDs{order(end-4:end)})

iM = 3;
iSubj = 5;
mnFR_nonZero = matFR_dT(iSubj).matSU_mnFR;
mnFR_nonZero(matFR_dT(iSubj).matSU_mnFR==0) = realmin;
tempStd = cellfun(@std, matFR_dT(iSubj).matSU_matFR, num2cell(zeros(size(matFR_dT(iSubj).matSU_matFR))), ...
    num2cell(ones(size(matFR_dT(iSubj).matSU_matFR)).*2), 'UniformOutput', false);
tempCV = cellfun(@rdivide, tempStd, mat2cell(mnFR_nonZero, [3000 3000 3000]', ones(1, size(tempStd, 2))), 'UniformOutput', false);

matIndLowCV=[];
for iCell = 1:length(matFR_dT(iSubj).setCellIDs)
    curCV = tempCV{iM, iCell};
    [y, ind] = sort(curCV, 'ascend');
    numZero = sum(curCV==0);
    indLowCV = ind(numZero+1:numZero+3000*0.02);
    pct = 0.02;
    matIndLowCV(:,iCell) = indLowCV;
    
    figure(100); clf;
    set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
    h(1) = subplot(5,1,1);
    imagesc(matFR_dT(iSubj).matSU_matFR{iM, iCell}');
    colormap(flipud(hot))
    set(h(1), 'CLim', [0 1])
    title(sprintf('%s: Cell %s, Movie %d', matFR_dT(iSubj).nameSubj, matFR_dT(iSubj).setCellIDs{iCell}, iM))
    h(2) = subplot(5,1,2);
    line([ind(numZero+1:numZero+3000*pct) ind(numZero+1:numZero+3000*pct)]', repmat([0 max(matFR_dT(iSubj).matSU_mnFR(3000*(iM-1)+1:3000*iM, iCell))]', 1, 3000*pct), 'Color', [1 1 1].*0.8)
    hold on;
    plot(h(2), matFR_dT(iSubj).matSU_mnFR(3000*(iM-1)+1:3000*iM, iCell), 'k-', 'LineWidth', 2);
    title('Averaged response across trials')
    h(3) = subplot(5,1,3);
    plot(h(3), tempStd{iM, iCell})
    title('Standard deviation across trials')
    h(4) = subplot(5,1,4);
    plot(h(4), tempCV{iM, iCell}, '.')
    title('Coefficient of variation across trials')
    h(5) = subplot(5,1,5);
    plot(h(5), 1./tempCV{iM, iCell}, 'r.')
    title('SNR across trials')
    set(h, 'XTick', 500:500:3000, 'XTickLabel', 50:50:300)
    set(h(2:5), 'TickDir', 'out', 'Box', 'off')
    input('')
end

% get the movie period where the cv was low
indImage = zeros(length(matFR_dT(iSubj).setCellIDs), 3000);
for iCell = 1:length(matFR_dT(iSubj).setCellIDs)
indImage(iCell, matIndLowCV(:,iCell)) = 1;
end

grandMnFR = mean(matFR_dT(iSubj).matSU_mnFR(3000*(iM-1)+1:3000*iM, :), 2);

figure;
subplot(2,1,1)
plot(sum(indImage, 1))
subplot(2,1,2)
plot(grandMnFR)
hold on;

indMov = find(sum(indImage,1)>2)';
plot(indMov, max(grandMnFR)+0.3, 'r*')

aaa = 1./cell2mat(tempCV);
iM = 3;
figure
imagesc(aaa(3000*(iM-1)+1:3000*iM, :)')
figure
hist(nansum(aaa(3000*(iM-1)+1:3000*iM,:), 2))
figure
plot(aaa(3000*(iM-1)+1:3000*iM,:), '-')
hold on
plot(nansum(aaa(3000*(iM-1)+1:3000*iM,:), 2), 'ko-')
% find(sum(aaa(3000*(iM-1)+1:3000*iM,:), 2)>50)
indMov = find(nansum(aaa(3000*(iM-1)+1:3000*iM,:), 2)>23);
cat(2, floor(indMov./600), mod(indMov./10, 60))



%% Across-neuron variation (correlation) change in time
iM = 1;
iSubj = 2;

tWindowStart = [0:0.5:299].*10;
tWindowEnd = [1:0.5:300]*10;

curMatFR_mn = matFR_dT(iSubj).matSU_mnFR(3000*(iM-1)+1:3000*iM, :);
curMatFR_norm = zscore(curMatFR_mn);
% 
% writerObj = VideoWriter(sprintf('%s_movie%d_betweenCellCorr2.avi', matFR_dT(iSubj).nameSubj, iM));
% open(writerObj);

numBin = length(tWindowStart);
matRho = NaN(length(matFR_dT(iSubj).setCellIDs), length(matFR_dT(iSubj).setCellIDs), length(tWindowStart));
for iBin = 1:numBin
    figure(2)
    [tempRho] = corr(curMatFR_norm(tWindowStart(iBin)+1:tWindowEnd(iBin), :), 'type', 'Spearman');
    imagesc(tempRho)
    set(gca, 'CLim', [-1 1])
    tempTime = cat(2, floor([tWindowStart(iBin)+1 tWindowEnd(iBin)]'./600), mod([tWindowStart(iBin)+1 tWindowEnd(iBin)]'./10, 60));
    title(sprintf('%s: Movie %d: %02d:%2.1f - %02d:%2.1f (bin %d/%d)', matFR_dT(iSubj).nameSubj, iM,...
        tempTime(1,1), tempTime(1,2), tempTime(2,1), tempTime(2,2), iBin, numBin))
%     input('')
    
    F = getframe;
    tempMovie(iBin) = F;
    matRho(:,:,iBin) = tempRho;
%     writeVideo(writerObj, F);
end

outputObj = VideoWriter(sprintf('/projects/parksh/NeuralBOLD/_labNote/_figs/%s_movie%d_betweenCellCorr_2hz.avi', ...
    matFR_dT(iSubj).nameSubj, iM), profile);
outputObj.FrameRate = 2;
open(outputObj)
for iFrame = 1:length(tempMovie)
writeVideo(outputObj, tempMovie(iFrame));
end
close(outputObj)
% close(writerObj);

%% Subcortical maps
flagBiowulf = 0;
if flagBiowulf
    dirDataHome = '/data/parks20/procdata/NeuroMRI/';
    addpath('/data/parks20/analysis/NeuroMRI/'); % to use doConv.m function
else
    ss = pwd;
    if ~isempty(strfind(ss, 'Volume')) % if it's local
        dirProjects = '/Volumes/PROJECTS';
        dirProcdata = '/Volumes/PROCDATA';
        dirDataHome = fullfile(dirProcdata, 'parksh');
        dirLibrary = '/Volumes/LIBRARY';
        addpath(fullfile(dirLibrary, 'matlab_utils'));
    else % on virtual machine
        dirProjects = '/projects';
        dirProcdata = '/procdata';
        dirDataHome = fullfile(dirProcdata, 'parksh');
        dirLibrary = '/library';
        addpath(fullfile(dirLibrary, 'matlab_utils'));
    end
end
setNameSubjNeural = {'Ava', 'Dav', 'Spi', 'Tor', 'Mat', 'Dan'}; % {'Tor', 'Rho', 'Sig', 'Spi'};
setNameArea = {'ML', 'ML', 'AF','AF', 'AM', 'AM'};
setMovie = [1 2 3]; %[4 5 6]; %[1 2 3];
tempS = num2str(setMovie);
MovieStr = tempS(~isspace(tempS));

% Get all the 
nameSubjBOLD = 'Art'; %'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);
load(fullfile(dirDataBOLD, sprintf('%sD99_ROIsegmentation.mat', nameSubjBOLD)), 'voltc_seg')

indSubcorticalArea(1).name = 'PI';
indSubcorticalArea(1).indROI = 196;
indSubcorticalArea(2).name = 'PL';
indSubcorticalArea(2).indROI = 197;
indSubcorticalArea(3).name = 'PM';
indSubcorticalArea(3).indROI = 198;
indSubcorticalArea(4).name = 'Bi';
indSubcorticalArea(4).indROI = 42;
indSubcorticalArea(5).name = 'Bmc';
indSubcorticalArea(5).indROI = 147;
indSubcorticalArea(6).name = 'Bpc';
indSubcorticalArea(6).indROI = 168;
indSubcorticalArea(7).name = 'Ld';
indSubcorticalArea(7).indROI = 9;
indSubcorticalArea(8).name = 'Lv';
indSubcorticalArea(8).indROI = 21;
indSubcorticalArea(9).name = 'LGNm';
indSubcorticalArea(9).indROI = 200;
indSubcorticalArea(10).name = 'LGNp';
indSubcorticalArea(10).indROI = 201;
indSubcorticalArea(11).name = 'Claustrum';
indSubcorticalArea(11).indROI = 176;
indSubcorticalArea(12).name = 'Striatum';
indSubcorticalArea(12).indROI = 175;


for iArea = 1:length(indSubcorticalArea)
    locArea = find(voltc_seg==indSubcorticalArea(iArea).indROI);
    indSubcorticalArea(iArea).locVox = locArea;
end

for iArea = 1:11 %length(indSubcorticalArea)
    matCorr = [];
    for iSubj = 2:6
        nameSubjNeural = setNameSubjNeural{iSubj};
        dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
        
        load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie%s_new_D99resample.mat', nameSubjNeural, nameSubjBOLD, MovieStr)),...
            'catVoltc_new')
        
        nCell = size(catVoltc_new,4);
        
        
        for iCell=1:nCell
            curVoltc = catVoltc_new(:,:,:,iCell);
            matCorr = cat(2, matCorr, curVoltc(indSubcorticalArea(iArea).locVox));
        end
        
    end
    
    % add pulvinar neurons
    nameSubjNeural = 'Dex';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
        
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie1_new_D99resample.mat', nameSubjNeural, nameSubjBOLD)),...
        'catVoltc_new')
    
    nCell = size(catVoltc_new,4);
        
    for iCell=1:nCell
        curVoltc = catVoltc_new(:,:,:,iCell);
        matCorr = cat(2, matCorr, curVoltc(indSubcorticalArea(iArea).locVox));
    end
    
    %
    catCorr{iArea} = matCorr;
    catCorr_mean{iArea} = mean(matCorr);
end

matSCCorr_mean = cell2mat(cat(1, catCorr_mean(:)));
[catAreaName{1:length(indSubcorticalArea)}] = deal(indSubcorticalArea.name);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(matSCCorr_mean)
set(gca, 'CLim', [-1 1].*0.4)
colorbar
ylabel('ROI')
xlabel('Cells')
nCell_area = [12 40 48 30 12 21];
set(gca, 'XTick', cumsum(nCell_area))
set(gca, 'YTick', 1:length(indSubcorticalArea), 'YTickLabel', catAreaName)
% print(gcf, fullfile(dirFig, 'PulvinarAmygdala_DavTorSpiMatDan_Movie123'), '-r150', '-dtiff')

%% Fingerprinting results
dirDataHome = '/procdata/parksh';
setNameSubjNeural = {'Ava', 'Dav', 'Spi', 'Mat', 'Dan'}; % {'Tor', 'Rho', 'Sig', 'Spi'};
setMovie = [1 2 3]; %[4 5 6]; %[1 2 3];
tempS = num2str(setMovie);
MovieStr = tempS(~isspace(tempS));

iSubj = 5; %5; %2; %5;

nameSubjNeural = setNameSubjNeural{iSubj};
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
switch lower(nameSubjNeural)
    case 'spi'
        nameSession_FPrint = 'Spice180120_other';
    case 'dav'
        nameSession_FPrint = 'Davida180515_other';
    case 'dan'
        nameSession_FPrint = 'Dango180123_other';
end
dirFPrint = fullfile(dirDataNeural, '_orgData', nameSession_FPrint, 'FPrint');

% get the file names of finger printing data
clear filename_fp cellID_fp d_fp
d_fp = dir(fullfile(dirFPrint, '*.mat'));
[filename_fp{1:length(d_fp)}] = deal(d_fp.name);
for iF = 1:length(filename_fp)
cellID_fp{iF} = filename_fp{iF}(1:strfind(filename_fp{iF}, '_')-1);
end

% get the list of selected cells
load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%sArtMovie%s_new.mat', nameSubjNeural, MovieStr)), 'paramCorr')

% going through relevant cells and get the finger printing results
matFR_fp = NaN(10, 6, size(paramCorr.validChanID, 1));
clear catDate
for iCell = 1:size(paramCorr.validChanID, 1)
    tt = strcmp(sprintf('%d', str2double(paramCorr.validChanID(iCell,:))), cellID_fp);
    indCell = find(tt);
    [curMat nDate] = response_calc_SHP(fullfile(dirFPrint, filename_fp{indCell}));
    matFR_fp(:, :, iCell) = curMat;
    catDate(iCell, 1) = nDate;
end
clear catMat
for iCell = 1:size(matFR_fp, 3)
    curMat = matFR_fp(:,:,iCell);
    curMat_norm = curMat./max(max(curMat));
    catMat(:,iCell) = curMat_norm(:);
end
figure
imagesc(catMat')
colormap(hot);
colorbar;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
title(sprintf('%s neurons: Fingerprinting normalized responses', nameSubjNeural))
ylabel('Cells')
xlabel('Stimulus')
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
print(gcf, fullfile(dirFig, sprintf('%s_fingerprinting_Movie%s', nameSubjNeural, MovieStr)), '-r150', '-dtiff')

matFR_fp_mean = squeeze(mean(matFR_fp));
figure; bar(matFR_fp_mean')

matFR_fp_faceother=[];
matFR_fp_faceother(1,:) = mean(matFR_fp_mean(1:2, :));
matFR_fp_faceother(2,:) = mean(matFR_fp_mean(3:6, :));
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1400 300])
bar(matFR_fp_faceother')
set(gca, 'XTick', 1:length(matFR_fp_faceother))
xlabel('Cell Index')
title(sprintf('%s: fingerprinting results: movie %s: faces vs. others', nameSubjNeural, MovieStr));
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
print(gcf, fullfile(dirFig, sprintf('%s_fingerprinting_Movie%s_avgFacesOthers', nameSubjNeural, MovieStr)), '-depsc')

% for movie 123
% for Spice: excChanInd = [6 11 30 35]; excChanID = {'09, '14', '42', '47'}  % Cells that should be excluded
% for Davida: indvalid = [1 4 5 8 9 10]; validchanID = paramCorr.validChanID(indvalid,:); %{'02' '08' '09' '17' '18' '20'}
% for Dango:  indvalid = [3 4 5 6 7 8 9];  validchanID =paramCorr.validChanID(indvalid,:) %{'10' '13' '15' '17' '18' '23' '26'}

% for movie 456
% for Spice: excChanInd = [4 7 8 9 14 15 18 23 32 34 38 41 42 43]; indvalid = setdiff(1:49, excChanInd);
% for Davida: indvalid = [2 3 4 9 10 12 14]; validchanID = paramCorr.validChanID(indvalid,:) %  {'02', '03', '07', '17', '18', 20', '22'};
% for Dango: indvalid = [1 4 5 6 7 8]; validchanID = paramCorr.validChanID(indvalid,:) % {'04', '09', '10', '17', '23', '26'};
%%
filename = '/procdata/parksh/Dan/matSDF_DavSpiMatDan_Movie123.mat'; % '/procdata/parksh/Dan/matSDF_AvaDavSpiMatDan_Movie456.mat'
dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';

load(filename)

%% Clustering

matFR_TR = cat(2, matSDF.matFR_SU_norm);
matNeuralRGR = cat(2, matSDF.matFR_SU);

numRepeat = 2; %10; %100; % number of repetition for entire clustering
[a, c, totalSS] = kmeans(matNeuralRGR', 1);
ClusteringSDF.totalSS = totalSS;
[a, c, totalSS] = kmeans(matFR_TR', 1);
ClusteringSDFnorm.totalSS = totalSS;

setK = 2:10;
opts = statset('Display','final');
numReplicates = 5; %30; %5;

for iK = 1:length(setK)
    K = setK(iK);
    % 3. Fine temporal resolution
    SU_indCluster_FR_SU = NaN(size(matNeuralRGR, 2), numRepeat);
    SU_sumD_FR_SU = NaN(K, numRepeat);
    % 3. Fine temporal resolution: Normalized
    SU_indCluster_FR_SU_norm = NaN(size(matFR_TR, 2), numRepeat);
    SU_sumD_FR_SU_norm = NaN(K, numRepeat);
    for iRep = 1:numRepeat
        fprintf(1, ':: K = %d; SDF ::\n', K);
        % 3. Fine temporal resolution
        [IDX_SUFR, C, SUMD_SUFR] = kmeans(matNeuralRGR', K, 'Replicates', numReplicates, 'Options', opts);
        SU_indCluster_FR_SU(:, iRep) = IDX_SUFR;
        SU_sumD_FR_SU(:, iRep) = SUMD_SUFR;
        % 4. Fine temporal resolution: Normalized
        [IDX_SUFRnorm, C, SUMD_SUFRnorm] = kmeans(matFR_TR', K, 'Replicates', numReplicates, 'Options', opts);
        SU_indCluster_FR_SU_norm(:, iRep) = IDX_SUFRnorm;
        SU_sumD_FR_SU_norm(:, iRep) = SUMD_SUFRnorm;
    end
    ClusteringSDF.resultKMeans(iK).SU_indCluster = SU_indCluster_FR_SU;
    ClusteringSDF.resultKMeans(iK).SU_sumD = SU_sumD_FR_SU;
    ClusteringSDFnorm.resultKMeans(iK).SU_indCluster = SU_indCluster_FR_SU_norm;
    ClusteringSDFnorm.resultKMeans(iK).SU_sumD = SU_sumD_FR_SU_norm;
end
figure
catSumD_norm = [];
for iK = 1:length(setK)
    catSumD_norm = cat(1, catSumD_norm, mean(mean(ClusteringSDFnorm.resultKMeans(iK).SU_sumD, 2)));
end
plot(catSumD_norm, 'o-')

matWSS=[];
matExpVar=[];
for iK = 1:length(setK)
curK = setK(iK);
matWSS(:,iK) = mean(mean(ClusteringSDFnorm.resultKMeans(iK).SU_sumD, 2));
end
totalSS = ClusteringSDFnorm.totalSS;
betweenSS = totalSS-matWSS;
propExplained = (totalSS-matWSS)./totalSS; %
fig3b=figure;
set(fig3b, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 410 290])
plot(setK, propExplained', 'k-'); hold on


iK=4; %2;
curInd = ClusteringSDFnorm.resultKMeans(iK).SU_indCluster(:,1);
[y, ind] = sort(curInd);
% catCellID=[];
% for iSubj=1:5
%     catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s_%s', matSDF(iSubj).area, matSDF(iSubj).nameSubj), length(matSDF(iSubj).setCellIDs), 1)), matSDF(iSubj).setCellIDs));
% end
% [y, ind] = sort(curInd);
tempCat = cat(2, num2str(y), char(catCellID(ind,:)));
figure
imagesc(matFR_TR(1:3000,ind)') %imagesc(matFR_SU_norm(:,ind)')
set(gca, 'CLim', [0 10])
colormap(hot)
set(gca, 'YTick', 1:1:size(matFR_TR, 2), 'YTickLabel', tempCat)
% title('QuickKMeans: 3Clusters (normalized 10hz SDF-AvaDavSpiMatDan-Movie456) ')
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% print(gcf, fullfile(dirFig, 'QuickKMeans_3Clusters_normSDF_AvaDavSpiMatDan_Movie456'), '-r150', '-dtiff')
% Hierarchical clustering
Z = linkage(matFR_TR','ward','euclidean');
c = cluster(Z, 'maxclust', 3);
[y, ind] = sort(c);
imagesc(matFR_TR(:,ind)')
set(gca, 'YTick', 1:1:size(matFR_TR, 2), 'YTickLabel', tempCat)
colormap(hot)
set(gca, 'CLim', [0 10])
tempCat = cat(2, num2str(y), char(catCellID(ind,:)));
set(gca, 'YTick', 1:1:size(matFR_TR, 2), 'YTickLabel', tempCat)
title('QuickHierarchicalClustering: 3Clusters (normalized 10hz SDF-AvaDavSpiMatDan-Movie456) ')
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% print(gcf, fullfile(dirFig, 'QuickHierarchicalClustering_3Clusters_normSDF_AvaDavSpiMatDan_Movie456'), '-r150', '-dtiff')

%% PCA of timeseries
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

matFR_TR =[];
matNeuralRGR=[];
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
    
    matFR_TR = cat(2, matFR_TR, matSDF(iSubj).matFR_SU_norm(:, indValid));
    matNeuralRGR = cat(2, matNeuralRGR, matSDF(iSubj).matFR_SU(:, indValid));
    
    catIndArea = cat(1, catIndArea, ones(length(indValid), 1).*setIndArea(iSubj));
    catIndSubj = cat(1, catIndSubj, ones(length(indValid), 1).*iSubj);
    catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', matSDF(iSubj).area, matSDF(iSubj).nameSubj), length(indValid), 1)), matSDF(iSubj).setCellIDs(indValid)));
        
    setSubjID = {'Ava', 'Dav', 'Spi', 'Mat', 'Dan'};
end


flagMerge = 1;
if iMovieSet == 1 && flagMerge
    load('/procdata/parksh/Spi/2016Nov_movie/matSDF_TorRhoSigSpiMovie123.mat') %% Merge earlier AF data
    % need to check later about Spice cells here (whether I need to exclude
    % some cells here or not)
    
    % Some info of the earlier dataset
    setSubj_2016 = {'Tor', 'Rho', 'Sig', 'Spi'};
    setIndSubj_2016 = [6 7 8 9];
    for iS = 1:length(matSDF)
        nCell_2016(iS) = length(matSDF(iS).setCellIDs);
        
        catIndArea = cat(1, catIndArea, ones(length(matSDF(iS).setCellIDs), 1)*2); % 2 for AF (1 for ML, 3 for AM)
        catIndSubj = cat(1, catIndSubj, ones(length(matSDF(iS).setCellIDs), 1).*setIndSubj_2016(iS));
        catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', 'AF', setSubj_2016{iS}), length(matSDF(iS).setCellIDs), 1)), matSDF(iS).setCellIDs));
    end
    
    matFR_TR = cat(2, matFR_TR, cat(2, matSDF.matFR_SU_norm));    
    matNeuralRGR = cat(2, matNeuralRGR, cat(2, matSDF.matFR_SU));
    
    setSubjID = cat(2, setSubjID, setSubj_2016);
end

    
% matFR_SU_norm = cat(2, matSDF.matFR_SU_norm);
% matFR_SU = cat(2, matSDF.matFR_SU);
% 
% setIndArea = [1 1 2 3 3]; %[0 1 2 3 3]; %for matSDF_DavSpiMatDan_Movie123.mat % [1 1 2 3 3]; %for matSDF_AvaDavSpiMatDan_Movie456.mat % 1 for ML, 2 for AF, 3 for AM
% catIndArea=[];
% catCellID=[];
% for iSubj = 1:5
%     catIndArea = cat(1, catIndArea, ones(length(matSDF(iSubj).setCellIDs), 1).*setIndArea(iSubj));
%     catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', matSDF(iSubj).area, matSDF(iSubj).nameSubj), length(matSDF(iSubj).setCellIDs), 1)), matSDF(iSubj).setCellIDs));
% end

[coeff_n, score_n, latent_n, tsquared_n, explained_n] = pca(matFR_TR', 'Economy', 'off', 'Centered', 'off');
cMap = [1 0 0; 0 1 0; 0 0 1].*0.8; % assign color (R, G, B) for each face patch

marker = {'', 'o', '^', 'v', '<', '>', 'square', 'diamond', '+', '*'};
setSize = fliplr(25:8:90);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [500 500 600 550])
scatter3(score_n(:,1), score_n(:,2), score_n(:,3), setSize(catIndSubj), cMap(catIndArea,:), 'LineWidth', 2);
% scatter3(score_n(:,1), score_n(:,2), score_n(:,3), 45, cMap(catIndArea,:), 'fill');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% text(score_n(:,1), score_n(:,2), score_n(:,3), catCellID);
% Rotate current 3d plot horizontally left and right
% % Example 3D plot
% figure
% [X,Y] = meshgrid(-3:.125:3);
% Z = peaks(X,Y);
% meshc(Z)
clear F
frames = 360; % sets the length of video
F(frames) = struct('cdata',[],'colormap',[]);

for h=1:frames
    view(h, 16);
%     camorbit(sin(h/360*2*pi),0)     % This makes the movie loopable. Divide by 10 to make shift smaller. Also could use cos(h/360*2*pi)
    %camorbit(0,sin(h/360*2*pi)/5) % Vertical shift up and down
    set(gca, 'XLim', [-60 80], 'YLim', [-100 100], 'ZLim', [-60 60])
    set(gca, 'XTick', -60:20:80, 'YTick', -100:20:100, 'ZTick', -60:20:60)
    xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
    drawnow
    F(h) = getframe(gcf);
end
% Replay the collected movie just for fun
fig = figure;
movie(fig,F,2)

% Write to video
v = VideoWriter(fullfile(dirFig, '3DPCA_10hzNormTS_Movie123_DavSpiMatDanTorRhoSigSpi_RotatingPlot.avi'));
v.Quality = 100; 
open(v);
writeVideo(v, F);
close(v);

D_SU_cosine = pdist(matFR_TR', 'cosine');
Z_SU_cosine = squareform(D_SU_cosine);
figure
imagesc(Z_SU_cosine)
colormap(flipud(hot));
set(gca, 'YTick', 1:length(catCellID), 'YTickLabel', catCellID)
set(gca, 'CLim', [0 1])
set(gca, 'CLim', [0 0.3])


%% PCA of correlation maps
% Set directories 
setNameSubjNeural = {'Ava', 'Dav', 'Spi', 'Mat', 'Dan'}; %{'Dav', 'Tor', 'Rho', 'Sig', 'Spi', 'Mat', 'Dan'}; %{'Tor', 'Rho', 'Sig', 'Spi'};
setNameArea = {'ML', 'ML', 'AF', 'AM', 'AM'};
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = '/procdata/parksh'; %'/data/parks20/procdata/NeuroMRI/'; %fullfile(dirProcdata, 'parksh');
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

setMovie = [4 5 6]; %[1 2 3]; %[4 5 6]; %[1 2 3]; %[4 5 6]; %[1 2 3];
tempS = num2str(setMovie);
MovieStr = tempS(~isspace(tempS));

% load the masks (movie-driven & brain-only mask)
load(fullfile(dirDataBOLD, sprintf('%s_MaskArrays.mat', nameSubjBOLD)), 'movieDrivenAmp'); %, 'brainMask_BlockAna3D');


% Load the corr map and select valid channel: subject by subject
numSubject = size(setNameSubjNeural, 2);
setIndArea = [1 1 2 3 3]; %[0 1 2 3 3]; %for matSDF_DavSpiMatDan_Movie123.mat % [1 1 2 3 3]; %for matSDF_AvaDavSpiMatDan_Movie456.mat % 1 for ML, 2 for AF, 3 for AM
catIndArea=[];
catCellID=[];
matR_SU_all = [];
for iSubj = 1:numSubject %2:5 %1:numSubject
    nameSubjNeural = setNameSubjNeural{iSubj}; %'Spi'; %'Tor';
    dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
    load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie%s_new.mat', nameSubjNeural, nameSubjBOLD, MovieStr)), 'matR_SU', 'paramCorr');
    
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
    
    matR_SU_valid = matR_SU(:, indValid);
    matR_SU_all = cat(2, matR_SU_all, matR_SU_valid);
    
    clear matR_SU matR_SU_valid
    
    catIndArea = cat(1, catIndArea, ones(length(indValid), 1).*setIndArea(iSubj));
    catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', setNameArea{iSubj}, nameSubjNeural), length(indValid), 1)), paramCorr.validChanID(indValid,:)));

end



[nx ny nz] = size(movieDrivenAmp.mask_amp1);
nVox = nx*ny*nz;

% Apply movie-driven mask to correlation matrix 
moviemask_vec = reshape(movieDrivenAmp.mask_amp1, nVox, 1); % change the 3D mask to 1D
matR_SU_all_moviemask = matR_SU_all(moviemask_vec,:); %matR_SU(moviemask_vec,:); % 15495 voxels

[coeff_map, score_map, latent_map, tsquared_map, explained_map] = pca(matR_SU_all_moviemask', 'Economy', 'off', 'Centered', 'off');
cMap = [1 0 0; 0 1 0; 0 0 1].*0.8; % assign color (R, G, B) for each face patch

figure;
scatter3(score_map(:,1), score_map(:,2), score_map(:,3), 45, cMap(catIndArea,:), 'fill');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
text(score_map(:,1), score_map(:,2), score_map(:,3), catCellID);

D_cosine = pdist(matR_SU_all_moviemask', 'cosine');
Z_cosine = squareform(D_cosine);
figure
imagesc(Z_cosine)
colormap(flipud(hot));
set(gca, 'YTick', 1:length(catCellID), 'YTickLabel', catCellID)
set(gca, 'CLim', [0 1])
set(gca, 'CLim', [0 0.3])



%% Focus on TR-resolution time series
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

matFR_TR =[];
matNeuralRGR=[];
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
    
    matFR_TR = cat(2, matFR_TR, matSDF(iSubj).matFR_TR(:, indValid));
    matNeuralRGR = cat(2, matNeuralRGR, matSDF(iSubj).matNeuralRGR(:, indValid));
    
    catIndArea = cat(1, catIndArea, ones(length(indValid), 1).*setIndArea(iSubj));
    catIndSubj = cat(1, catIndSubj, ones(length(indValid), 1).*iSubj);
    catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', matSDF(iSubj).area, matSDF(iSubj).nameSubj), length(indValid), 1)), matSDF(iSubj).setCellIDs(indValid)));
        
    setSubjID = {'Ava', 'Dav', 'Spi', 'Mat', 'Dan'};
end


flagMerge = 1; %0; %1;
if iMovieSet == 1 && flagMerge
    load('/procdata/parksh/Spi/2016Nov_movie/matSDF_TorRhoSigSpiMovie123.mat') %% Merge earlier AF data
    % need to check later about Spice cells here (whether I need to exclude
    % some cells here or not)
    
    % Some info of the earlier dataset
    setSubj_2016 = {'Tor', 'Rho', 'Sig', 'Spi'};
    setIndSubj_2016 = [6 7 8 9];
    for iS = 1:length(matSDF)
        nCell_2016(iS) = length(matSDF(iS).setCellIDs);
        
        catIndArea = cat(1, catIndArea, ones(length(matSDF(iS).setCellIDs), 1)*2); % 2 for AF (1 for ML, 3 for AM)
        catIndSubj = cat(1, catIndSubj, ones(length(matSDF(iS).setCellIDs), 1).*setIndSubj_2016(iS));
        catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', 'AF', setSubj_2016{iS}), length(matSDF(iS).setCellIDs), 1)), matSDF(iS).setCellIDs));
    end
    
    matFR_TR = cat(2, matFR_TR, cat(2, matSDF.matFR_TR));    
    matNeuralRGR = cat(2, matNeuralRGR, cat(2, matSDF.matNeuralRGR));
    
    setSubjID = cat(2, setSubjID, setSubj_2016);
end

indValid = ~isnan(matNeuralRGR);
% matNeuralRGR_noNaN = matNeuralRGR(indValid);
matNeuralRGR_noNaN = reshape(matNeuralRGR(indValid), 375-21, size(matNeuralRGR, 2));

D = pdist(matNeuralRGR_noNaN'); %, 'cosine');
Z = squareform(D);
figure
imagesc(Z)
colormap(flipud(hot));
set(gca, 'YTick', 1:length(catCellID), 'YTickLabel', catCellID)

[matR_neuralRGR] = corr(matNeuralRGR, 'rows', 'complete', 'type', 'Spearman');
figure
imagesc(matR_neuralRGR)
set(gca, 'CLim', [-1 1])
set(gca, 'YTick', 1:length(catCellID), 'YTickLabel', catCellID)
[matR_TR] = corr(matFR_TR, 'rows', 'complete', 'type', 'Spearman');
figure
subplot(1,2,1)
imagesc(matR_neuralRGR)
set(gca, 'CLim', [-1 1])
subplot(1,2,2)
imagesc(matR_TR)
set(gca, 'CLim', [-1 1])
title('TR-resolution')
axis square
set(gca, 'YTick', 1:length(catCellID), 'YTickLabel', catCellID)
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
% print(gcf, fullfile(dirFig, 'pairwiseSpearmanCorr_2p4sec_movie123_cLim-1to1'), '-dtiff', '-r150')

%%% Looking into pair-wise correlations in detail
% [matR_SU] = corr(matFR_SU_norm, 'rows', 'complete', 'type', 'Spearman');
aa = tril(matR_TR, -1);
[i, j] = find(aa > 0.8);
cat(2, catCellID(i), catCellID(j))
aa_vector = aa(:);
[sortedR, sortedSub] = sort(aa_vector, 'descend');
figure
plot(sortedR, 'o-')
sortedR(1:10)
sortedR(10:30)
length(find(sortedR>0.6))
length(find(sortedR>0.7))
help hist
bin  = -1:0.1:1;
[n, indBin] = histc(aa_vector, bin);
find(indBin==17)
[ii, jj] = ind2sub(size(aa), find(indBin==17))
cat(2, catCellID(ii), catCellID(jj))
curBin = 16;
[ii, jj] = ind2sub(size(aa), find(indBin==curBin))
cat(2, catCellID(ii), catCellID(jj))
strfind(catCellID, 'AFTOr069a')
tloc = strfind(catCellID, 'AFTor069a')
find(cellfun('isempty', strfind(catCellID, 'AFTor069a'))<1)
plot(aa(:, 78), 'o-')
[sortR78, indsortcell78] = sort(matR_SU(78,:), 'descend');
sortR78(1:10)
indsortcell78(1:10)
catCellID(indsortcell78(1:10))

% PCA on the 2.4s time series
[coeff_n, score_n, latent_n, tsquared_n, explained_n] = pca(zscore(matFR_TR', 'Economy', 'off', 'Centered', 'off');
cMap = [1 0 0; 0 1 0; 0 0 1].*0.8; % assign color (R, G, B) for each face patch

marker = {'', 'o', '^', 'v', '<', '>', 'square', 'diamond', '+', '*'};
setSize = fliplr(25:8:90);

figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [500 500 600 550])
scatter3(score_n(:,1), score_n(:,2), score_n(:,3), setSize(catIndSubj), cMap(catIndArea,:), 'LineWidth', 2);
% scatter3(score_n(:,1), score_n(:,2), score_n(:,3), 45, cMap(catIndArea,:), 'fill');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
text(score_map(:,1), score_map(:,2), score_map(:,3), catCellID);


%% 2019/05/22 after saving the TS from multiple subject for each area
dirDataHome = '/procdata/parksh/_macaque/';
setMovie = [1 2 3];
setArea = {'AF', 'AM'};

% get all the time courses across areas
matFR_TR = []; matNeuralRGR=[];
catCellID=[];
catIndArea=[];
catIndSubj=[];
for iArea = 1:length(setArea)
    nameArea = setArea{iArea}; %'AF';
    load(fullfile(dirDataHome, sprintf('matSDF_%s_Movie%s.mat', nameArea, num2str(setMovie, '%d%d%d'))))
    
    clear tempCurCat
    tempCurCat = cat(2, matSDF.matFR_TR);
    matFR_TR = cat(2, matFR_TR, tempCurCat); %cat(2, matFR_TR, matSDF(iSubj).matFR_TR(:, indValid));
    matNeuralRGR = cat(2, matNeuralRGR, cat(2, matSDF.matNeuralRGR));
    
    catIndArea = cat(1, catIndArea, ones(size(tempCurCat, 2), 1).*iArea);
    %     catIndSubj = cat(1, catIndSubj, ones(length(indValid), 1).*iSubj);
    for iSubj = 1:length(matSDF)
        catIndSubj = cat(1, catIndSubj, ones(length(matSDF(iSubj).setCellIDs), 1).*iSubj + 10*iArea);
        catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', matSDF(iSubj).area, matSDF(iSubj).nameSubj), length(matSDF(iSubj).setCellIDs), 1)), matSDF(iSubj).setCellIDs));
    end
end

flagMerge = 1; %0; %1;
if flagMerge %iMovieSet == 1 && flagMerge
    load('/procdata/parksh/_macaque/Spi/2016Nov_movie/matSDF_TorRhoSigSpiMovie123.mat') %% Merge earlier AF data
    % need to check later about Spice cells here (whether I need to exclude
    % some cells here or not)
    
    % Some info of the earlier dataset
    setSubj_2016 = {'Tor', 'Rho', 'Sig', 'Spi'};
    setIndSubj_2016 = [3 4 5 6];
    for iS = 1:length(matSDF)
        nCell_2016(iS) = length(matSDF(iS).setCellIDs);
        
        catIndArea = cat(1, catIndArea, ones(length(matSDF(iS).setCellIDs), 1)*1); % 1 for AF (2 for AM)
        catIndSubj = cat(1, catIndSubj, ones(length(matSDF(iS).setCellIDs), 1).*setIndSubj_2016(iS) + 10*1);
        catCellID = cat(1, catCellID, strcat(cellstr(repmat(sprintf('%s%s', 'AF', setSubj_2016{iS}), length(matSDF(iS).setCellIDs), 1)), matSDF(iS).setCellIDs));
    end
    
    matFR_TR = cat(2, matFR_TR, cat(2, matSDF.matFR_TR));    
    matNeuralRGR = cat(2, matNeuralRGR, cat(2, matSDF.matNeuralRGR));
    
%     setSubjID = cat(2, setSubjID, setSubj_2016);
end
% indValid = ~isnan(matNeuralRGR);
% % matNeuralRGR_noNaN = matNeuralRGR(indValid);
% matNeuralRGR_noNaN = reshape(matNeuralRGR(indValid), 375-21, size(matNeuralRGR, 2));


[matR_TR] = corr(matFR_TR, 'rows', 'complete', 'type', 'Spearman');
figure
imagesc(matR_TR)
set(gca, 'CLim', [-1 1])
title('TR-resolution')
axis square
set(gca, 'YTick', 1:length(catCellID), 'YTickLabel', catCellID)

% PCA on the 2.4s time series
[coeff_n, score_n, latent_n, tsquared_n, explained_n] = pca(zscore(matFR_TR)', 'Economy', 'off', 'Centered', 'on');
cMap = [1 0 0; 0 1 0; 0 0 1].*0.8; % assign color (R, G, B) for each face patch

marker = {'', 'o', '^', 'v', '<', '>', 'square', 'diamond', '+', '*'};
setSize = fliplr(25:8:90);
subjectCMap{1} = flipud(jet(20));
subjectCMap{2} = jet(10);

figPCA= figure;
set(figPCA, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [500 500 600 550]);
% scatter3(score_n(:,1), score_n(:,2), score_n(:,3), setSize(catIndSubj), cMap(catIndArea,:), 'LineWidth', 2);
scatter3(score_n(:,1), score_n(:,2), score_n(:,3), 45, cMap(catIndArea,:), 'fill');
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
text(score_n(:,1), score_n(:,2), score_n(:,3), num2str(catIndSubj));


figPCA= figure;
set(figPCA, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [500 500 1100 950]);
for iArea = 1:2
    indValidArea = find(catIndArea==iArea);
    figure(figPCA);
%     plot3(score_n(indValidArea,1), score_n(indValidArea,2), score_n(indValidArea,3), 'o', 'MarkerSize', 15,...
%         'MarkerFaceColor', subjectCMap{iArea}(catIndSubj(indValidArea) - 10*iArea, :), 'MarkerEdgeColor', subjectCMap{iArea}(catIndSubj(indValidArea) - 10*iArea, :));
%     hold on;
    scatter3(score_n(indValidArea,1), score_n(indValidArea,2), score_n(indValidArea,3), 45,...
        subjectCMap{iArea}(catIndSubj(indValidArea) - 10*iArea, :), 'fill');
    hold on;
end
figure(figPCA);
text(score_n(:,1), score_n(:,2), score_n(:,3), num2str(catIndSubj));
axis tight
xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
view([0 90])
% scatter3(score_n(:,1), score_n(:,2), score_n(:,3), 45, subjectCMap(catIndSubj,:), 'fill');

dim = [1 2];
figPCA2D=figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position',  [46    56   928   919])
for iArea = 1:2
    indValidArea = find(catIndArea==iArea);
    figure(figPCA2D);
    scatter(score_n(indValidArea, dim(1)), score_n(indValidArea, dim(2)), 50, subjectCMap{iArea}(catIndSubj(indValidArea) - 10*iArea, :),...
        'fill') %'LineWidth', 2)
    hold on;
end
text(score_n(:,dim(1)), score_n(:,dim(2)), num2str(catIndSubj));
axis tight
xlabel(sprintf('PC%d', dim(1))); 
ylabel(sprintf('PC%d', dim(2))); 


% compare PCs between cells and fMRI
 load('/procdata/parksh/_macaque/Art/Art_movieTS_fMRI_Movie123_PCA.mat', 'resultsPCA', 'pc1_fMRI');
 
TR=2.4;
k = gampdf([-40:TR:40],4,2);
addpath('/library/matlab_utils');


convCoeff=[];
for iMov = 1:length(indMovieNeuron)
    curCoeff = coeff_n((iMov-1)*125+8:iMov*125, 1:5); %S(validC(iChan), indMovieNeuron(iMov)).mnFR
    curCoeff = curCoeff-repmat(mean(curCoeff), 118, 1); % centering
    curCoeff = doConv(curCoeff,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
    curCoeff = cat(1, NaN(7,5), curCoeff'); %curNeuralTC(1:7) = NaN;
    convCoeff = cat(1, convCoeff, curCoeff); % concatenation across movies
end




% % Select voxels: voxels that showed correlation higher than 0.3 with any one
% % of the neurons
% critCorr = 0.3; 
% critNumNeuron = 6; %13;
% matValidVox = abs(matR_SU_all_moviemask)>critCorr;
% locValidVox = sum(matValidVox, 2)>critNumNeuron;
% matR_SU_all_moviemask_valid  = matR_SU_all_moviemask(locValidVox, :);