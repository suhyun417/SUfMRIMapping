% genFig_SFN2015
%
% make a variety of figures for SFN 2015 poster

%% Params
% color-blind friendly colors
color_sb = [86 180 233]./255; % sky blue
color_bg = [0 158 115]./255; % bluish green
color_o = [230 159 0]./255; % orange
color_rp = [204 121 167]./255; % reddish purple

dirFig = '/projects/parksh/NeuralBOLD/_labNote/_figs/';
% dirFig = '/Volumes/PROJECTS/parksh/NeuralBOLD/_labNote/_figs/';
% 
addpath('/library/matlab_utils/')


%% Sample cell SDF to show diversity
dirDataNeural = '/procdata/parksh/Tor'; %'/Volumes/PROCDATA/parksh/Tor';
S = createCellRegressor_indMov(dirDataNeural, {'065a', '082a', '089a'}, 1, 1000); % last input is sigma (for SDF computation) in msec

tempC = bone(15); % jet(size(S(1).matsdf{1}, 2))*0.75; % temporary ColorOrder for line plot

figSDF = figure;
set(figSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 600])
% case of two example neurons
SP(1) = subplot('position', [0.1 0.55 0.8 0.35]);
SP(2) = subplot('position', [0.1 0.15 0.8 0.35]);
set(SP, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(SP(1), S(1).matsdf{1}, 'LineWidth', 2)
plot(SP(2), S(3).matsdf{1}, 'LineWidth', 2)
set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
set(SP, 'YLim', [0 30], 'YTick', [0 30], 'YTickLabel', [0 30], 'XTickLabel', [])
% save
print(figSDF, fullfile(dirFig, 'exampleSDF_cell082a089a_trialbytrial'), '-depsc')

% For DAL's BSC review: 5 cells
dirDataNeural = '/procdata/parksh/Tor'; %'/Volumes/PROCDATA/parksh/Tor';
setCellIDs =  {'065a', '075a', '082a', '089a', '118a'};
S = createCellRegressor_indMov(dirDataNeural, {'065a', '075a', '082a', '089a', '118a'}, 1, 1000); % last input is sigma (for SDF computation) in msec

tempC = bone(15); % jet(size(S(1).matsdf{1}, 2))*0.75; % temporary ColorOrder for line plot

figSDF = figure;
set(figSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 800 860])
% case of five example neurons
startY = 0.07;
height = 0.15;
gap = 0.03;
for iCell = 1:5
    SP(iCell) = subplot('position', [0.1 startY+height*(iCell-1)+gap*(iCell-1) 0.8 height]);
end
set(SP, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
for iCell = 1:5
    plot(SP(iCell), S(iCell).matsdf{1}, 'LineWidth', 2)
end

set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
set(SP(3:5), 'YLim', [0 30], 'YTick', [0 30], 'YTickLabel', [0 30])
set(SP(1:2), 'YLim', [0 10], 'YTick', [0 10], 'YTickLabel', [0 10])
set(SP(2:5), 'XTickLabel', [])
set(SP(1), 'XTickLabel', 0:50:300)
xlabel(SP(1), 'Time (sec)')
ylabel(SP(3), 'Spike / sec')
% save
print(figSDF, fullfile(dirFig, sprintf('exampleSDF_cell%s_trialbytrial', cell2mat(setCellIDs))), '-depsc')

% Legend
figLegend = figure;
set(figLegend, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 200 50])
set(gca, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
x = 1:11; a = 0.5;
line([x; x+a], [zeros(size(x)); zeros(size(x))+tan(pi/8)], 'LineWidth', 5)
xlim([0.5 11.5]);
print(figLegend, fullfile(dirFig, 'exampleSDF_cell082a089a_trialbytrial_figLegend'), '-depsc')


%% Methods
% Sample neural regressor
setMovie = [1 2 3];
indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);

k = gampdf([-40:2.4:40],4,2); % MION function

neuralrgrs=[]; 
iChan = 11;
neuralrgrs = cat(1, S(iChan, indMovieNeuron).mnFR);
neuralrgrs = neuralrgrs-mean(neuralrgrs); % centering
neuralrgrs = doConv(neuralrgrs,k); % convolve MION kernel %conv(neuralrgrs,k,'same');

% Sample voxel TS
voxCoords(1,:) = [7 33 19];
voxCoords(2,:) = [7 18 10];
% voxCoords(3,:) = [7    28    25];

dirDataBOLD = '/procdata/parksh/Art';
load(fullfile(dirDataBOLD, 'Art_movieTS_fMRI_indMov.mat'))

% Get movie IDs common in two dataset
commonSetMovie = intersect(paramBOLD.unimov, paramSDF.setMovIDs);
dataBOLD.mvoltc = voltcIndMov(commonSetMovie);
dataBOLD.unimov = commonSetMovie;

setMovie = [1 2 3];

% 1. fMRI tc
fmritc=[];
indMovieBOLD = find(ismember(dataBOLD.unimov, setMovie)>0);
for iM = indMovieBOLD %1:length(indMovieBOLD)
    curvoltc = dataBOLD.mvoltc{iM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc = cat(4,fmritc,pcvoltc);
end
clear dataBOLD voltcIndMov

% Compare voxel TS and neural regressor
figExCorr = figure;
set(figExCorr, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 1200 600])
for ivox=1:size(voxCoords, 1)

SP(ivox) = subplot(size(voxCoords,1), 1, ivox);
tempTS = squeeze(fmritc(voxCoords(ivox,1), voxCoords(ivox,2), voxCoords(ivox,3),:)).*(-1); % MION
zscoreTempTS =  (tempTS-nanmean(tempTS))./nanstd(tempTS);
plot(zscoreTempTS, 'k-', 'LineWidth', 3);
hold on
plot(zscore(neuralrgrs), 'm-', 'LineWidth', 3)
axis tight

end
set(SP, 'YLim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [], 'XTickLabel', [])
set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
print(figExCorr, fullfile(dirFig, 'exampleVoxTS_cell082a_0p58_0p28'), '-depsc')
    

%% Map of averaged SDF
dirDataNeural = '/procdata/parksh/Tor';
load(fullfile(dirDataNeural, 'Tor_movieTS_SU_indMov.mat'))

setMovie = [1 2 3];
validC = [1:8, 11:34, 36:51]; % valid channel with movie [1 2 3]
indMovieNeuron = find(ismember(paramSDF.setMovIDs, setMovie)>0);
% MION function
TR=2.4;
k = gampdf([-40:TR:40],4,2);

for iChan = 1:length(validC) % compute correlation channel-by-channel
    tempFR = [];
    tempFR = cat(1, S(validC(iChan), indMovieNeuron).mnFR);
    matFR_TR(:,iChan) = tempFR;
end
meanFR_TR = mean(matFR_TR, 2); % averaged of all single units

% make a regressor
neuralrgrs=meanFR_TR;
neuralrgrs = neuralrgrs-mean(neuralrgrs); % centering
neuralrgrs = doConv(neuralrgrs,k); % convolve MION

% Compute a correlation map
[nx, ny, nz, nt] = size(fmritc);
nVox = nx*ny*nz;
 [Rvals, Pvals] = corr(reshape(fmritc, nVox, nt)', neuralrgrs',...
        'rows','complete', 'type', 'Spearman');    
mapR_averagedSDF = reshape(Rvals, [nx, ny, nz]).*(-1); % because of MION







% % Reshape BOLD 4-d data
% [nx, ny, nz, nt] = size(fmritc);
% nVox = nx*ny*nz;
% 
% matR_SU = NaN(nVox, length(validC)); matP_SU = NaN(nVox, length(validC));

% % individual neurons 
% dirDataNeural = '/procdata/parksh/Tor'; %'/Volumes/PROCDATA/parksh/Tor';
% FR_dT = createCellRegressor_indMov_discreteTime(dirDataNeural, {'065a', '082a', '089a', '128a'}, [1 2 3], 2.4);
% 
% figSDF = figure;
% set(figSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 700 1200 300])
% plot(S(1,1)
% % case of two example neurons
% % SP(1) = subplot('position', [0.1 0.55 0.8 0.35]);
% % SP(2) = subplot('position', [0.1 0.15 0.8 0.35]);
% % set(SP, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
% plot(SP(1), S(2).matsdf{1}, 'LineWidth', 2)
% plot(SP(2), S(3).matsdf{1}, 'LineWidth', 2)
% set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
% set(SP, 'YLim', [0 30], 'YTick', [0 30], 'YTickLabel', [0 30], 'XTickLabel', [])
% % save
% print(figSDF, fullfile(dirFig, 'exampleSDF_cell082a089a_trialbytrial'), '-depsc')
% 
% % %averaged across trials
% % dirDataNeural = '/Volumes/PROCDATA/parksh/Tor';
% % S = createCellRegressor_indMov(dirDataNeural, {'065a'}, 1);
% % meanSDF = S(1).mnsdf;
% % steSDF = std(S(1).matsdf{1}, 0, 2)./sqrt(size(S(1).matsdf{1},2)-1);
% % xAxis = 1:1:length(meanSDF);
% % figAvgSDF = figure;
% % set(figAvgSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 300])
% % % subplot('Position', [0 0 1 1]); cla;
% % p = patch([xAxis fliplr(xAxis)]', [meanSDF+steSDF; flipud(meanSDF-steSDF)], [1 1 1].*0.7);
% % set(p, 'EdgeColor', 'none'); 
% % hold on;
% % plot(S(1).mnsdf, 'k-', 'LineWidth', 4)
% % set(gca, 'XTick', [], 'YTick', [], 'Box', 'off', 'LineWidth', 3)
% % % save
% % print(figAvgSDF, fullfile(dirFig, 'exampleTS_cell065a_avgSDF'), '-depsc')


% % Sample resampled & convolved SDF
% k = gampdf([-40:2.4:40], 4, 2);
% resampleSDF = resample(meanSDF, 0.001*1000, 2.4*1000);
% convSDF = conv(resampleSDF, k, 'same');
% 
% figConvSDF = figure;
% set(figConvSDF, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 300])
% % subplot('Position', [0 0 1 1]); cla;
% plot(convSDF, 'ko', 'LineWidth', 4, 'MarkerSize', 12, 'MarkerFaceColor', 'w')
% xlim([0 125])
% set(gca, 'XTick', [], 'YTick', [], 'Box', 'off', 'LineWidth', 3)
% % save
% print(figConvSDF, fullfile(dirFig, 'exampleTS_cell065a_convSDF'), '-depsc')
% 
% 
% % Sample convolved principal component
% load('/Volumes/PROCDATA/parksh/Tor/eigen/pcaCatMovie123_SDF_tor.mat')
% pc1 = pcaMovCat.coeff(:,1);
% resamplePC1 = resample(pc1, 0.1*10, 2.4*10);
% convPC1 = conv(resamplePC1, k, 'same');
% 
% figConvPC1 = figure;
% set(figConvPC1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 920 175])
% plot(convPC1, 'ko', 'LineWidth', 4, 'MarkerSize', 6, 'MarkerFaceColor', 'w')
% xlim([0 375])
% set(gca, 'XTick', [], 'YTick', [], 'Box', 'off', 'LineWidth', 2)
% % save
% print(figConvPC1, fullfile(dirFig, 'exampleTS_PC1_convPC'), '-depsc')


%%
% Sample Voxel TS
voxCoords(1,:) = [30    10    19];
voxCoords(2,:) = [34 34 21];
voxCoords(3,:) = [7    28    25];

load('/Volumes/PROCDATA/parksh/dataNeuralBOLD_TorArt_indMov.mat', 'dataBOLD')

fmritc = cat(4, dataBOLD.catmvoltc{1:3});
for ivox=1:3
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1000 200])
% subplot('Position', [0 0 1 1])
plot(squeeze(fmritc(voxCoords(ivox,1), voxCoords(ivox,2), voxCoords(ivox,3),:)),...
    'Color', 'k', 'LineWidth', 4)
axis tight
set(gca, 'LineWidth', 2, 'box', 'off')
set(gca, 'YTick', [],'XTick', [])
print(gcf, sprintf(fullfile(dirFig, 'exampleTS_vox%d'), ivox), '-depsc')
end
% print(gcf, sprintf(fullfile(dirFig, 'movieRGR_%s'), char(varnames(id))), '-depsc'


%% Sample cell rasters
sampleCell = {'065a', '075a', '118a', '082a'}; %'/procdata/parksh/Tor/
win = [0 300]; % seconds
iMov = 1;

figRaster = figure;
set(figRaster, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 200])
for iCell = 1:length(sampleCell)
    fName = sprintf('/Volumes/PROCDATA/parksh/Tor/tmov%dsig%s.mat', iMov, char(sampleCell{iCell}));
    load(fName)
    
    switch lower(dat.h.units)
        case 'sec'
            win = [0 300]; % seconds or milliseconds
        case 'ms'
            win = [0 300*1000];
    end
    spikes = {};
    for t=1:length( dat.s)
        ts = dat.s{t};
        spikes{t} = ts(find((ts>=win(1)) & ts<=win(2)));
    end
    
    figure(figRaster); clf;
    subplot('Position', [0 0 1 1])
    iS=1;
    set(gca,'YLim',[iS-1 iS]);
        line([spikes{iS} spikes{iS}]', repmat([iS-1 iS]', 1, length(spikes{iS})), 'Color', 'k')
        
    for iS=1:length(spikes)
        set(gca,'YLim',[iS-1 iS]);
        line([spikes{iS} spikes{iS}]', repmat([iS-1 iS]', 1, length(spikes{iS})), 'Color', 'k')
        hold on;
        %     yline(spikes{iS},prop,val);
    end
    
    set(gca,'YLim',[0 length(dat.s)]);
    set(gca,'XLim', win);
    set(gca,'YTick', [], 'XTick', [])
    % save as eps
    print(figRaster, sprintf(fullfile(dirFig, 'exampleTS_cell%s'), char(sampleCell{iCell})), '-depsc')
    
end




%%
% Different movie regressors
setMovID = [1 2 3];
flagSM = 1; % flag for compression and smoothing

fullRGR4fps = createMovieRGR_4fps_indMov(setMovID, flagSM); %createFullMovieRegressors_4fps_indMov(setMovID); %
indValidRGR = [1, 3, 9, 20, 21, 22, 25]; 
% 1: 'Luminance', 3: 'Speed', 9: 'Beta Contrast', 20: 'Faces', 21: 'One face', 22: 'Body parts', 25: 'Any animal'
varnames = cat(1, fullRGR4fps(1).features(indValidRGR), {'Face size'});
indSample = [1 2 4]; % 8];

% Face scale regressor (in TR unit)
ttt=load('/Volumes/PROCDATA/parksh/MovieRegressors/dbtmMriReg.mat');
scaleRGR = ttt.reg.xx(7,:)';

catRGR=[];matRGR=[];
for iMov=1:length(setMovID)
    m = setMovID(iMov);
    matCurRGR = fullRGR4fps(m).smoRegressors(:,indValidRGR); %fullRGR4fps(iMov).regressors(:,indValidRGR); %fullRGR4fps(iMov).regressors;
    catRGR = cat(1, catRGR, matCurRGR); % concatenation across movies
end
matRGR = resample(catRGR, 0.25*100, 2.4*100);
matRGR = cat(2, matRGR, scaleRGR);

for iR = 1:length(indSample)
    id = indSample(iR);
    figure;
    set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 200])
%     subplot('Position', [0 0 1 1])
    plot(catRGR(:,id), 'Color', color_sb, 'LineWidth', 4)
    axis tight
    set(gca, 'LineWidth', 2, 'box', 'off')
    set(gca, 'YTick', [],'XTick', 1:400:3600, 'XTickLabel', [], 'TickDir', 'out')
    print(gcf, sprintf(fullfile(dirFig, 'movieRGR_%s'), char(varnames(id))), '-depsc')
end
id = 8; % face scale regressor
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1200 200])
% subplot('Position', [0 0 1 1])
plot(1:2.4:900, matRGR(:,id), 'Color', color_sb, 'LineWidth', 4)
axis tight
set(gca, 'LineWidth', 2, 'box', 'off')
set(gca, 'YTick', [],'XTick', 0:100:900, 'XTickLabel', [],'TickDir', 'out')
print(gcf, sprintf(fullfile(dirFig, 'movieRGR_%s'), char(varnames(id))), '-depsc')


%% PC and movie regressors
% PC1 and face size
scale_norm = (matRGR(:,8)-nanmean(matRGR(:,8)))./nanstd(matRGR(:,8));
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [110         307        1036         278])
plot(zscore(pcaMovCat.coeff(:,1)), 'k-', 'LineWidth', 2)
hold on;
plot(scale_norm, 'Color', color_bg, 'LineWidth', 2)
set(gca, 'LineWidth', 2, 'box', 'off', 'FontSize', 12, 'FontName', 'Arial')
axis tight
set(gca, 'XTick', 0:125:375, 'XTickLabel', [],'YTickLabel', [], 'TickDir', 'out')
print(gcf, fullfile(dirFig, 'pca_movieRGR_PC1_facesize'), '-depsc')

% PC3 and speed
% scale_norm = (matRGR(:,8)-nanmean(matRGR(:,8)))./nanstd(matRGR(:,8));
figure;
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [110         307        1036         278])
plot(zscore(pcaMovCat.coeff(:,3)), 'k-', 'LineWidth', 2)
hold on;
plot(zscore(matRGR(:,2)), 'Color', color_rp, 'LineWidth', 2)
set(gca, 'LineWidth', 2, 'box', 'off', 'FontSize', 12, 'FontName', 'Arial')
axis tight
set(gca, 'XTick', 0:125:375, 'XTickLabel', [],'YTickLabel', [], 'TickDir', 'out')
print(gcf, fullfile(dirFig, 'pca_movieRGR_PC3_speed'), '-depsc')

% 
figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 278 278])
plot(zscore(pcaMovCat.coeff(:,1)), scale_norm, 'o', 'MarkerEdgeColor', color_bg, 'MarkerFaceColor', 'w',...
    'MarkerSize', 10, 'LineWidth', 2)
set(gca, 'LineWidth', 2, 'box', 'off')
set(gca, 'XTickLabel', [],'YTickLabel', [], 'TickDir', 'out')
axis square
print(gcf, fullfile(dirFig, 'pca_movieRGR_PC1_facesize_scatter'), '-depsc')

figure
set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 278 278])
plot(zscore(pcaMovCat.coeff(:,3)), zscore(matRGR(:,2)), 'o',...
    'MarkerSize', 10, 'LineWidth', 2, 'MarkerEdgeColor', color_rp, 'MarkerFaceColor', 'w')
set(gca, 'LineWidth', 2, 'box', 'off')
set(gca, 'XTickLabel', [],'YTickLabel', [], 'TickDir', 'out')
axis square
print(gcf, fullfile(dirFig, 'pca_movieRGR_PC3_speed_scatter'), '-depsc')


