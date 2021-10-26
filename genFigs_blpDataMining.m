% genFigs_blpDataMining.m

flagLocal = 0; %1; % 1 for desktop

switch flagLocal
    case 1 % desktop
        dirData = '/Volumes/PROCDATA/parksh/';
        dirProject = '/Volumes/PROJECTS/parksh/NeuralBOLD/analysis/';
    case 0
        dirData = '/procdata/parksh/';
        dirProject = '/projects/parksh/NeuralBOLD/analysis/';
end

dirDataNeural = fullfile(dirData, 'Tor');
% load(fullfile(dirDataNeural, 'BLPRGR_indMov_tor.mat'))

% compute BLP in different temporal resolutions (to compare 
% differently sampled movie regressors)
% BLPRGR_TR = createCellRegressor_BLP_indMov_discreteTime(dirDataNeural, 2.4, 0);
% matMeanBLP_TR = cat(1,BLPRGR_TR.meanBLP);
% catBLP_TR=[];
% for iBP=1:4
%     catBLP_TR(:,iBP) = cat(1, matMeanBLP_TR{:,iBP});
% end
% 
% BLPRGR_4fps = createCellRegressor_BLP_indMov_discreteTime(dirDataNeural, 0.25);
% matMeanBLP_4fps = cat(1,BLPRGR_4fps.meanBLP);
% catBLP_4fps=[];
% for iBP=1:4
%     catBLP_4fps(:,iBP) = cat(1, matMeanBLP_4fps{:,iBP});
% end

BLPRGR_10fps = createCellRegressor_BLP_indMov_discreteTime(dirDataNeural, 0.1);
matMeanBLP_10fps  = cat(1,BLPRGR_10fps .meanBLP);
catBLP_10fps =[];
for iBP=1:4
    catBLP_10fps (:,iBP) = cat(1, matMeanBLP_10fps {:,iBP});
end

%% Compare BLP with SU activities
% 1. compare to averaged SDF of all SUs
% list of channel IDs and movie IDs of each file
d_n = dir(fullfile(dirDataNeural, '*sig*.mat'));
listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')'; % list of channels

setCellIDs = unique(listSU_all);
setMovIDs = [1 2 3];
sizeTimeBin_sec = 1/10; % just to match it to BLP in 10 fps
FR_dT = createCellRegressor_indMov_discreteTime(dirDataNeural, setCellIDs, setMovIDs, sizeTimeBin_sec);
FR_dT_sec = createCellRegressor_indMov_discreteTime(dirDataNeural, setCellIDs, setMovIDs, 1);

catSU = []; catSU_sec = [];
for iCell =1:size(FR_dT,1)
    if isempty(FR_dT(iCell,1).mnFR)
       continue;
    else
        tempCat = []; tempCat_sec = [];
        tempCat = cat(1, FR_dT(iCell,:).mnFR); % concatenate movies
        catSU = cat(2, catSU, tempCat); % Matrix of timeseries of all cells (time x cell)
        
        tempCat_sec = cat(1, FR_dT_sec(iCell,:).mnFR); % concatenate movies
        catSU_sec = cat(2, catSU_sec, tempCat_sec); % Matrix of timeseries of all cells (time x cell)
    end
end
meanSpikePerSec_SU = mean(catSU_sec); % mean firing rate (spikes/sec) for each SU
cellID=cellstr(cat(1,FR_dT(:,1).cellID));

% 1-1. compare to all SUs
% compute simple correlations, as always
R_indSU = corr(catBLP_10fps, catSU, 'rows', 'complete', 'type', 'Spearman');

% PLOT
figCorrBLPSU = figure;
set(figCorrBLPSU, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(R_indSU');
set(gca, 'YTick', 1:size(R_indSU,2), 'YTickLabel', cellID);
set(gca, 'XTick', 1:size(R_indSU,1), 'XTickLabel', 
set(gca, 'FontSize', 15)
% set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', fullRGR4fps(1).features(indValidRGR));
xlabel('Band-limited powers')
ylabel('Cell ID')

% make blue-white-red colorbar
cval = 0.6;
cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end

endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));

figure(figCorrBLPSU)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;

% save as graphic file
dirFig =  [dirProject, '../_labNote/_figs/'];
print(gcf, fullfile(dirFig, 'BLP_catMovie_corrIndSU'), '-dtiff', '-r150')

% 1-2. compare to first-order statistics of population (mean and variance)
% compute simple correlations, as always
R_avgSU = corr(catBLP_10fps, avgSU, 'rows', 'complete', 'type', 'Spearman');
R_varSU = corr(catBLP_10fps, varSU, 'rows', 'complete', 'type', 'Spearman');

% PLOT
figCorrBLPSU_mean = figure;
set(figCorrBLPSU_mean, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(R_avgSU');
set(gca, 'XTick', 1:size(R_avgSU,1), 'XTickLabel', BLPRGR_10fps(1).blpfreq(:,1))
set(gca, 'FontSize', 15)
% set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', fullRGR4fps(1).features(indValidRGR));
xlabel('Band-limited powers')

figure(figCorrBLPSU_mean)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;

figCorrBLPSU_var = figure;
set(figCorrBLPSU_var, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(R_varSU');
set(gca, 'XTick', 1:size(R_varSU,1), 'XTickLabel', BLPRGR_10fps(1).blpfreq(:,1))
set(gca, 'FontSize', 15)
% set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', fullRGR4fps(1).features(indValidRGR));
xlabel('Band-limited powers')

figure(figCorrBLPSU_var)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;

% 1-3. High-gamma BLP and high-firing neurons
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto')
plot(R_indSU(4,:), meanSpikePerSec_SU, 'ko', 'LineWidth', 2)
axis square
set(gca, 'FontSize', 15)
xlabel('Correlation between high-gamma BLP and single unit response')
ylabel('Firing rate (spikes/s)')

%% Compare BLP with movie features

% Create movie regressors
setMovID = [1 2 3];

flagSM = 1; % flag for compression and smoothing
fullRGR4fps = createMovieRGR_4fps_indMov(setMovID, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
indValidRGR = [1, 3, 9, 20, 21, 22, 25]; % 1: 'Luminance', 3: 'Speed', 9: 'Beta Contrast', 20: 'Faces', 21: 'One face', 22: 'Body parts', 25: 'Any animal'

% Face scale regressor (in TR unit)
ttt=load([dirData 'MovieRegressors/dbtmMriReg.mat']);
scaleRGR = ttt.reg.xx(7,:)';

catMovieRGR=[];
for iMov=1:length(setMovID)
    m = setMovID(iMov);
    matCurRGR = fullRGR4fps(m).regressors; %fullRGR4fps(m).smoRegressors(:,indValidRGR); %fullRGR4fps(m).regressors(:,indValidRGR); %fullRGR4fps(m).regressors;
    catMovieRGR = cat(1, catMovieRGR, matCurRGR); % concatenation across movies
end

[R_blp_fullRGR] = corr(catMovieRGR, catBLP_4fps, 'rows', 'complete');

varnames_fullRGR = fullRGR4fps(1).features; %cat(1, fullRGR4fps(1).features, {'Face size'}); %cat(1, fullRGR4fps(1).features(indValidRGR), {'Face size'});

figCorrBLPmovie_catMov = figure;
set(figCorrBLPmovie_catMov, 'Color', 'w', 'PaperPositionMode', 'auto')
imagesc(R_blp_fullRGR);
set(gca, 'YTick', 1:size(catMovieRGR,2), 'YTickLabel', varnames_fullRGR);
set(gca, 'FontSize', 15)
% set(gca, 'YTick', 1:length(indValidRGR), 'YTickLabel', fullRGR4fps(1).features(indValidRGR));
xlabel('Band-limited powers')
ylabel('Features')

% make blue-white-red colorbar
cval = 0.5;
cmin = -cval; cmax = cval;
colornum = 256;
colorInput = [1 0 0; 1 1 1; 0 0 1];
oldSteps = linspace(-1, 1, length(colorInput));
newSteps = linspace(-1, 1, colornum);
for j=1:3 % RGB
    newmap_all(:,j) = min(max(transpose(interp1(oldSteps, colorInput(:,j), newSteps)), 0), 1); 
end

endPoint = round((cmax-cmin)/2/abs(cmin)*colornum);
newmap = squeeze(newmap_all(1:endPoint, :));

figure(figCorrBLPmovie_catMov)
set(gca, 'CLim', [cmin cmax])
colormap(flipud(newmap))
set(gca, 'TickDir', 'out')
box off
c=colorbar;

dirFig =  [dirProject, '../_labNote/_figs/'];
print(gcf, fullfile(dirFig, 'BLP_catMovie_movieRGR'), '-dtiff', '-r150')


% Make a soundtrack of BLP

BLPRGR_30fps = createCellRegressor_BLP_indMov_discreteTime(dirDataNeural, 1/30);
matMeanBLP_30fps  = cat(1,BLPRGR_30fps .meanBLP);
catBLP_30fps =[];
for iBP=1:4
    catBLP_30fps (:,iBP) = cat(1, matMeanBLP_30fps {:,iBP});
end

for iMov=1:3
for iBP=1:4
inputSig = BLPRGR_30fps(iMov).meanBLP{iBP};
Fs_org = 30; % in samples per sec
for flagAM = 0:1 %1; % amplitude modulation
if flagAM    
saveFileName = fullfile(dirDataNeural, 'BLP_ost', sprintf('BLP_0%d_movie%d_AM.wav', iBP, iMov));
else 
    saveFileName = fullfile(dirDataNeural, 'BLP_ost', sprintf('BLP_0%d_movie%d_FM.wav', iBP, iMov));
end
contSignal2Soundtrack(inputSig, Fs_org, flagAM, saveFileName);
end
end
end
 

