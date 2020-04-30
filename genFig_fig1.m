% genFig_fig1.m
%
% Generate figure 1, which is about the methods to make a single-unit
% correlation map
% Example unit time course for each step, example voxel time courses,
% correlation map example

%% Settings
ss = pwd;
if ~isempty(strfind(ss, 'Volume')) % if it's local
    dirProjects = '/Volumes/PROJECTS';
    dirProcdata = '/Volumes/PROCDATA';
    dirLibrary = '/Volumes/LIBRARY';
else % on virtual machine
    dirProjects = '/projects';
    dirProcdata = '/procdata';
    dirLibrary = '/library';
end
    
% Add necessary toolbox % Should be 
addpath(fullfile(dirLibrary, 'matlab_utils')) % for convolution
% addpath(fullfile(dirProjects, 'parksh/_toolbox/Boot_Time_Series'))

% Set directories 
nameSubjNeural = 'Tor';
nameSubjBOLD ='Art'; % 'Ava'; %'Art'; % 'Ava'; %'Art'; %'Ava'; %'Art';
dirDataHome = fullfile(dirProcdata, 'parksh');
dirDataNeural = fullfile(dirDataHome, nameSubjNeural);
dirDataBOLD = fullfile(dirDataHome, nameSubjBOLD);

% Directory for saving figures as graphic files
dirFig = fullfile(dirProjects, 'parksh/NeuralBOLD/_labNote/_figs/');


%% Load the original time series
filenameNeural = [nameSubjNeural, '_movieTS_SU_indMov.mat'];
% fileNameNeural_BLP = [nameSubjNeural, '_movieTS_BLPLFP_indMov.mat'];
filenameBOLD = [nameSubjBOLD, '_movieTS_fMRI_indMov.mat'];

fprintf(1, '\nLoading single unit data of %s: %s ....', nameSubjNeural, filenameNeural)
load(fullfile(dirDataNeural, filenameNeural))
% load(fullfile(dirDataNeural, fileNameNeural_BLP))
fprintf(1, '\nLoading fMRI data of %s: %s ....\n', nameSubjBOLD, filenameBOLD)
load(fullfile(dirDataBOLD, filenameBOLD))


%% Example neural time course (Fig 1a)
indUnit = find(strcmp('082a', paramSDF.setCellIDs)>0);
 
% 1. downsampled & averaged across trials
tempC = bone(15); % jet(size(S(1).matsdf{1}, 2))*0.75; % temporary ColorOrder for line plot

fig1a1 = figure;
set(fig1a1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1000 300])
set(gca, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(S(indUnit, 1).matFR{1}, 'o-', 'LineWidth', 2)
hold on
plot(S(indUnit, 1).mnFR, 'ko-', 'LineWidth', 4)
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XLim', [0 125], 'XTick', 0:25:125, 'XTickLabel', [])
set(gca, 'YLim', [0 70], 'YTick', [0 70], 'YTickLabel', [0 70])
% % save
% print(fig1a1, fullfile(dirFig, 'fig1a_1'), '-depsc')


% 2. HRF convolved
k = gampdf([-40:2.4:40],4,2); % MION function

neuralrgrs=[]; 
indUnit = find(strcmp('082a', paramSDF.setCellIDs)>0);
neuralrgrs = S(indUnit, 1).mnFR;
neuralrgrs = neuralrgrs-mean(neuralrgrs); % centering
neuralrgrs = doConv(neuralrgrs,k); % convolve MION kernel %conv(neuralrgrs,k,'same');

fig1a2 = figure;
set(fig1a2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1000 200])
% set(gca, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(neuralrgrs, 'ko-', 'LineWidth', 2); hold on;
% line([125 250; 125 250], [-6 -6; 6 6], 'LineStyle', '--', 'Color', 'k')
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XLim', [0 125], 'XTick', 0:25:125, 'XTickLabel', [])
set(gca, 'YLim', [-6 6], 'YTick', [-6 0 6], 'YTickLabel', [])
% save
print(fig1a2, fullfile(dirFig, 'fig1a_2'), '-depsc')

% 3. HRF convolved
k = gampdf([-40:2.4:40],4,2); % MION function

indUnit = find(strcmp('082a', paramSDF.setCellIDs)>0);
neuralrgrs=[];
for iMov = 1:3
    curNeuralTC = S(indUnit, iMov).mnFR;
    curNeuralTC = curNeuralTC-mean(curNeuralTC); % centering
    curNeuralTC = doConv(curNeuralTC,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
    curNeuralTC(1:7) = NaN;
    
    neuralrgrs = cat(2, neuralrgrs, curNeuralTC); % concatenation across movies
end
% neuralrgrs = cat(1, S(indUnit, [1 2 3]).mnFR);
% neuralrgrs = neuralrgrs-mean(neuralrgrs); % centering
% neuralrgrs = doConv(neuralrgrs,k); % convolve MION kernel %conv(neuralrgrs,k,'same');
% setOnsetTR = [1:7, 125:125+6, 250:250+6];
% neuralrgrs(setOnsetTR) = NaN;

fig1a3 = figure;
set(fig1a3, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1000 200])
% set(gca, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(neuralrgrs, 'ko-', 'LineWidth', 2); hold on;
line([125 250; 125 250], [-6 -6; 6 6], 'LineStyle', '--', 'Color', 'k')
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XLim', [0 375], 'XTick', 0:25:375, 'XTickLabel', [])
set(gca, 'YLim', [-6 6], 'YTick', [-6 0 6], 'YTickLabel', [])
% save
print(fig1a3, fullfile(dirFig, 'fig1a_3'), '-depsc')


%% Example voxel time course (Fig 1b)
% 1. raw time series from each trial & averaged across trials
% Needs to be loaded after launching BlockAna and run S_neuralRegressor
cd(fullfile(dirProjects, '/parksh/NeuralBOLD/analysis/BlockAna'))
blockana;
S_neuralRegressor;
% Subject name and directory
nameSubjBOLD = 'Art';
sessionFileList = 'FL_e66_allfiles2.txt'; % 'FL_ava_allfiles2.txt';
skip=7; % number of TRs to skip

SI = SU_createSessInfo(sessionFileList,[],skip); 
filelist = SU_makeFileList(SI.monkID,SI.sessID,SI.scanID);

MAX_TR = SI.max_tr;
skip=SI.skiptr;

unimov = [1 2 3];
nunimov = length(unimov);
catimgdat = [];
% rgr = zeros(1,nunimov*MAX_TR);
% procmovs = [];

for u = 1:nunimov
    movindxs = find(SI.movID == unimov(u));
    nmovindxs = length(movindxs);
    
    % collect all the data files for this movie
    filelist = SU_makeFileList(SI.monkID(movindxs),SI.sessID(movindxs),...
        SI.scanID(movindxs));
    
    totfiles = length(filelist);
    
    for f=1:totfiles
        s_sub = filelist{f}.subj;
        s_ses = filelist{f}.sess;
        s_sc  = filelist{f}.scan;
        
        [fmri_tc,dgz,notes] = SU_loadFile(s_sub,s_ses,s_sc);        
        
        if f==1
            [a,b,c,t] = size(fmri_tc);
            allvoltc = zeros(a,b,c,MAX_TR,totfiles);
        end
        
        scanlen = size(fmri_tc,4);
        if scanlen == 250
            mov_start_tr = 63;% offset when movie starts in scan
        elseif scanlen == 125
            mov_start_tr = 1; % short movie
        else
            fprintf('WARNING: bad scan length %d\n',scanlen);
        end
        
        val_tr = [mov_start_tr:(mov_start_tr+MAX_TR-1)];
        [clp_fmri_tc,clp_mdgz] = SU_clipMovDat(fmri_tc,dgz,val_tr);
        allvoltc(:,:,:,:,f) = clp_fmri_tc;
        alldgz(f) = dgz;
    end
    
end

voxCoords(1,:) = [7 33 19];
voxCoords(2,:) = [7 18 10];

% voxCoords(1,:) = [30    10    19];
% voxCoords(2,:) = [34 34 21];
% voxCoords(3,:) = [7    28    25];

voxtc_pc=[];
for iVox = 1:size(voxCoords,1)%2; %1;
    voxtc = squeeze(allvoltc(voxCoords(iVox,1), voxCoords(iVox,2), voxCoords(iVox,3), :, :)); 
    tempvoxtc_pc = ((voxtc-repmat(mean(voxtc), 125, 1))./repmat(mean(voxtc), 125, 1))*100;
    
    voxtc_pc = cat(3, voxtc_pc, -tempvoxtc_pc(8:125,:));
end

tempC = cool(60); %autumn(40); %cool(60); %spring(40); % jet(size(S(1).matsdf{1}, 2))*0.75; % temporary ColorOrder for line plot

%vox1
fig1b1 = figure;
set(fig1b1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1000 300])
set(gca, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(squeeze(voxtc_pc(:,:,1)), 'o-', 'LineWidth', 2); 
hold on
plot(mean(voxtc_pc(:,:,1),2), 'mo-', 'LineWidth', 4)
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XLim', [0 125]-7, 'XTick', [0:25:125]-7, 'XTickLabel', [])
set(gca, 'YLim', [-20 20], 'YTick', [-20:10:20], 'YTickLabel', [])
% save
print(fig1b1, fullfile(dirFig, 'fig1b_1_vox1'), '-depsc')

% vox2
fig1b1 = figure;
set(fig1b1, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1000 300])
set(gca, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(squeeze(voxtc_pc(:,:,2)), 'o-', 'LineWidth', 2); 
hold on
plot(mean(voxtc_pc(:,:,2),2), 'mo-', 'LineWidth', 4)
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XLim', [0 125]-7, 'XTick', [0:25:125]-7, 'XTickLabel', [])
set(gca, 'YLim', [-15 15],  'YTick', [-15:5:15], 'YTickLabel', [])
% save
print(fig1b1, fullfile(dirFig, 'fig1b_1_vox2'), '-depsc')

% Legend
figLegend = figure;
set(figLegend, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 200 50])
% subplot('Position', [0 0 1 1])
set(gca, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
x = 1:36; a = 1.5;
line([x; x+a], [zeros(size(x)); zeros(size(x))+tan(pi/8)], 'LineWidth', 2)
xlim([0 38]);
axis off
%save
print(figLegend, fullfile(dirFig, 'fig1b_1_legend'), '-depsc')

   
% 2. concatenated across movies
setMovie = [1 2 3];
fmritc=[];
for iM = setMovie %1:length(indMovieBOLD)
    curvoltc = voltcIndMov{iM};
    avgvoltc = repmat(nanmean(curvoltc,4),[1 1 1 size(curvoltc,4)]);
    if ~isempty(find(avgvoltc==0, 1))
        avgvoltc(avgvoltc==0) = realmin; % get rid of zeros because it causes NaNs in percent signals
    end
    pcvoltc = ((curvoltc - avgvoltc)./avgvoltc)*100;
    fmritc = cat(4,fmritc,pcvoltc);
end

% vox1
fig1b2 = figure;
set(fig1b2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1000 200])
% set(gca, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(-squeeze(fmritc(voxCoords(1,1),voxCoords(1,2),voxCoords(1,3),:)) , 'mo-', 'LineWidth', 2); hold on;
line([125 250; 125 250], [-1 -1; 1 1]*10, 'LineStyle', '--', 'Color', 'k')
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XLim', [0 375], 'XTick', [0:25:375], 'XTickLabel', [])
set(gca, 'YLim', [-1 1]*10, 'YTick', [-1 0 1]*10, 'YTickLabel', [])
% save
print(fig1b2, fullfile(dirFig, 'fig1b_2_vox1'), '-depsc')

% vox2
fig1b2 = figure;
set(fig1b2, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [100 100 1000 200])
% set(gca, 'ColorOrder', tempC, 'NextPlot', 'replacechildren')
plot(-squeeze(fmritc(voxCoords(2,1),voxCoords(2,2),voxCoords(2,3),:)) , 'mo-', 'LineWidth', 2); hold on;
line([125 250; 125 250], [-1 -1; 1 1]*10, 'LineStyle', '--', 'Color', 'k')
set(gca, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
set(gca, 'XLim', [0 375], 'XTick', [0:25:375], 'XTickLabel', [])
set(gca, 'YLim', [-1 1]*10, 'YTick', [-1 0 1]*10, 'YTickLabel', [])
% save
print(fig1b2, fullfile(dirFig, 'fig1b_2_vox2'), '-depsc')


%% Compute correlation: example single-unit map with sample voxel time series (Fig 1c)
%
voxCoords(1,:) = [7 33 19]; % rho = 0.58
voxCoords(2,:) = [7 18 10]; % rho = 0.28

fig1c = figure;
set(fig1c, 'Color', 'w', 'PaperPositionMode', 'auto', 'Position', [600 600 1200 600])
for ivox=1:size(voxCoords, 1)
    figure(fig1c);
    SP(ivox) = subplot(size(voxCoords,1), 1, ivox);
    tempTS = squeeze(fmritc(voxCoords(ivox,1), voxCoords(ivox,2), voxCoords(ivox,3),:)).*(-1); % MION
    zscoreTempTS =  (tempTS-nanmean(tempTS))./nanstd(tempTS);
    plot(zscoreTempTS, 'm-', 'LineWidth', 3);
    hold on
    plot(zscore(neuralrgrs), 'k-', 'LineWidth', 3)
    axis tight
end
set(SP, 'XLim', [0 375], 'XTick', 0:25:375, 'XTickLabel', [])
set(SP, 'YLim', [-3 3], 'YTick', [-3 0 3], 'YTickLabel', [], 'XTickLabel', [])
set(SP, 'TickDir', 'out', 'Box', 'off', 'LineWidth', 2, 'FontSize', 15)
% save
print(fig1c, fullfile(dirFig, 'fig1c_tc'), '-depsc') 


% Example sagittal slice with two voxels annotated
cd /projects/parksh/NeuralBOLD/analysis/BlockAna/
blockana;
S_neuralRegressor;

load(fullfile(dirDataNeural, sprintf('CorrMap_SU_%s%sMovie123.mat', nameSubjNeural, nameSubjBOLD)), 'paramCorr', 'matR_SU') 
indUnit = find(strcmp('082a', cellstr(paramCorr.validChanID))>0);

[nx, ny, nz, nt] = size(fmritc);
nVox = nx*ny*nz;
mapR082a_3d= reshape(matR_SU(:,indUnit), [nx, ny, nz]); % because of MION

global DSP
DSP.proc.scalarmap_3d = mapR082a_3d;



