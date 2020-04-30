function [S, dataBOLD] = getCorrMap(nameSubjNeural, nameSubjBOLD, saveFileName_matSDF, saveFileName)
% clear all;

%% Environmental setup
% global flagLocal

flagLocal = 0; %2; %0; %1; % 0 for server (Pauli), 1 for local machine, 2 for home

switch flagLocal
    case 0
    addpath('/einstein0/share/UCNI_Library/matlab_utils/');
    dirData = '/einstein0/USRlab/data/parksh/'; %/rmov/';
    case 1
    addpath('/Volumes/share/UCNI_Library/matlab_utils/');
    dirData = '/Volumes/USRlab/data/parksh/'; %/rmov/';
    case 2
    addpath('/Users/shyuny/Documents/MATLAB/expToolbox/BlockAna/'); %'/einstein0/share/UCNI_Library/matlab_utils/');
    addpath('/Users/shyuny/Documents/MATLAB/expToolbox/BlockAna/matlab_utils');
    dirData = '/Volumes/HONGKONG/USRlab/data/parksh/'; %'/Users/shyuny/ResearchProjects/0NeuralBOLDCorr/data/';
end
% W = what(dirData);

% Subject name and directory
nameSubjNeural = 'Tor'; %'Rho';
nameSubjBOLD = 'Art';
dirDataNeural = [dirData, nameSubjNeural, '/'];
dirDataBOLD = [dirData, nameSubjBOLD, '/'];


%% Get the list of Single Units and Movie IDs 
% Get the movie list from the fMRI data
d_b = dir(fullfile(dirDataBOLD, '*BOLD.mat')); %W = what(dirDataBOLD);
load(fullfile(dirDataBOLD, d_b(1).name)); % for now it's only one file
listMovBOLD = MCD.unimov; % Movie IDs of the current BOLD

% Get the movie list from single unit data
d_n = dir(fullfile(dirDataNeural, '*sig*.mat'));
[listMatSUFile{1:length(d_n), 1}] = deal(d_n.name); % Get the list of single unit files 
% listMatSUFile = regexp(cat(2, d_n.name), '\w*.mat', 'match')'; %regexp(cat(2, d_n.name), '(?<=mat)\w*', 'match')';
% listSUchannelID = unique(regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match'));
listMovSU = unique(regexp(cat(2,d_n.name), '\d*(?=sig)', 'match')); % Movie IDs of the single units in cellstr

% Choose the movies, or select common movie IDs between single units and
% BOLD as much as possible
setMovIDs = input('Movie IDs (e.g. [1 2 3])? (default: common movies in single units and fMRI)');
if isempty(setMovIDs)
    setMovIDs = intersect(listMovBOLD, str2num(char(listMovSU))');
end
% setMovIDs = [13 14 15]; %intersect(listMovBOLD, str2num(char(listMovSU))');
[~, locMovBOLD] = intersect(listMovBOLD, setMovIDs);
% [setMovIDs, locMovBOLD] = intersect(listMovBOLD, str2num(char(listMovSU))'); 

% Get the list of single units & their movie
listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')';
listSUchannelID = unique(listSU_all);

indMov=[];
for iCh=1:length(listSUchannelID) % go through channel-by-channel     
     tempListMov{iCh} = regexp([listMatSUFile{strcmp(listSUchannelID{iCh}, listSU_all)}], '\d*(?=sig)', 'match');
     % Get indices for common movies across cells
     indMov(iCh, :) = ismember(setMovIDs, str2num(char(tempListMov{iCh}))'); % data presence matrix (Channels x MovieIDs) 
end

% List of cells that have data for the movies
setCellIDs = listSUchannelID(sum(indMov,2)==length(setMovIDs));

% setCellIDs = {'006a'; '009a'; '010a'; '012a'; '013a'}; % Rhombus %{'003a';'005a';'007a';'009a';'013b';'014a'}; % Sig: '013a' has very low spikes
% setMovIDs  = [7 8 9 10 11 12]; % Rhom %[7 8 9]; % Sig % [1 2 3 7 8 9];


%% Load fMRI data
fprintf(1, '\nPreparing BOLD time series...\n');
% fprintf(1, 'Number of movies in original time series: %d\n', length(setMovIDs))

dataBOLD = MCD; % Get necessary params for further DSP in BlockAna: for now, just all the params

dataBOLD = rmfield(dataBOLD, 'catmvoltc'); % remove the previous time series

lengthIndMovie_TR = dataBOLD.max_tr;
% lengthIndMovie_sec = lengthIndMovie_TR*dataBOLD.TR; % 300;

% Get the necessary part of time courses, corresponding to our (common) movie set
% if ~issorted(setMovIDs) % just to be sure
%     setMovIDs = sort(setMovIDs);
% end

countTR = 0; tc=[];
for iMov = 1:length(locMovBOLD) % concatenate across different movies
    fprintf(1, 'Movie ID: %d (%d/%d), Location in original BOLD: %d\n', setMovIDs(iMov), iMov, length(setMovIDs), locMovBOLD(iMov));
    tc(:,:,:,countTR+1:countTR+lengthIndMovie_TR) = MCD.catmvoltc(:,:,:,lengthIndMovie_TR*(locMovBOLD(iMov)-1)+1:lengthIndMovie_TR*locMovBOLD(iMov)); % should be optimized for different set of movies
    countTR = size(tc, 4);
end
dataBOLD.catmvoltc = tc;
fprintf(1, 'Time course for %d movies has been saved in dataBOLD.tc\n', length(setMovIDs));
fprintf(1, 'The number of TRs of dataBOLD.catmvoltc should be %d and is %d\n', ...
    lengthIndMovie_TR*length(setMovIDs), size(dataBOLD.catmvoltc, 4));

dataBOLD.unimov = setMovIDs;

clear MCD tc %SigData





%% Get mean SDF and regressors for each single unit
S = createCellRegressor(dirDataNeural, setCellIDs, setMovIDs); 
% -- S: for each cell (S.cellID), 
%       1) matsdf: for each movie (cell), matSDF for each trial, 2) mnsdf: mean SDF across trials, 
%       3) catmnsdf: concatenated mean SDF across movies
%       4) rgrsMION: MION function convolved regressors (before resampling)
%       5) rgrsMION_resample: MION function convolved regressors after
%       resampling (in accordance with the number of frames in fMRI data)

% 
fig_mnsdf = figure; clf;
set(fig_mnsdf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Name', 'Mean SDF across trials')

fig_rgrsMIONresample = figure; clf;
set(fig_rgrsMIONresample, 'Color', 'w', 'PaperPositionMode', 'auto', 'Name', 'MION convolved regressors')

sp=[]; sp2=[];
for iC = 1:length(S)
  
    fprintf(1, 'Channel %s (%d/%d)', S(iC).cellID, iC, length(S));
    figure(fig_mnsdf);
    %   sp(iC) = subplot(length(S),1,iC);
    plot(S(iC).catmnsdf);
    axis tight;
    %   title(sprintf('cell: %s, movies: [%d %d %d]', S(iC).cellID, setMovIDs))
    
    figure(fig_rgrsMIONresample);
    %   sp2(iC) = subplot(length(S),1,iC);
    plot(S(iC).rgrsMION_resample);
    axis tight;
    %   title(sprintf('cell: %s, movies: [%d %d %d]', S(iC).cellID, setMovIDs))
    input('')
    
end
  
% set(sp, 'XTickLabel', [])
% set(sp(end), 'XTickLabel', 0:100:900)
% xlabel(sp(end), 'Time (s)')
% ylabel(sp(end), 'Firing rate (spikes/s)')
% 
% set(sp2, 'XTick', 0:25:375, 'XTickLabel', [], 'XLim', [0 375])
% set(sp2(end), 'XTickLabel', [0:25:375].*2.4)
% xlabel(sp2(end), 'Time (s)')
% ylabel(sp2(end), 'Activity (a.u.)')




%% Replace initial stimulus-onset responses with NaNs (in accordance with fMRI responses)
onsetResp_numTR = dataBOLD.skiptr; %7;
onsetResp_sec = onsetResp_numTR*dataBOLD.TR; %tr;
onsetResp_msec = onsetResp_sec*1000;

% 
figure;
set(fig_mnsdf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Name', 'MION regressors and BOLD response from one example voxel')
sp3=[];
for iC = 1:length(S)
    fprintf(1, 'Channel %s (%d/%d)', S(iC).cellID, iC, length(S));
    for iM = 1:length(setMovIDs)
        S(iC).rgrsMION_resample(lengthIndMovie_TR*(iM-1)+1:lengthIndMovie_TR*(iM-1)+onsetResp_numTR) = NaN;
    end
%     sp3(iC) = subplot(length(S)+1,1,iC);
    plot(S(iC).rgrsMION_resample);
%     input('')
%     title(sprintf('cell: %s, movies: [%d %d %d]', S(iC).cellID, setMovIDs))    
end
% sp3(iC+1) = subplot(length(S)+1, 1, length(S)+1);
% plot(1:375,squeeze(dataBOLD.tc(22,33,20,:)), 'r-')
% title('BOLD time course from one example voxel')
% 
% set(sp3, 'XTick', 0:25:375, 'XTickLabel', [], 'XLim', [0 375])
% set(sp3(end), 'XTickLabel', [0:25:375].*2.4)
% xlabel(sp3(end), 'Time (s)')
% ylabel(sp3(end), 'fMRI activity (% BOLD)')



%% Compute correlation
% Reshape BOLD 4-d data
[nx, ny, nz, nt] = size(dataBOLD.catmvoltc);
nVox = nx*ny*nz;

fprintf(1, 'Compute correlation using regressors:: \n')
for iC = 1:length(S)
    tic
    [Rvals, Pvals] = corr(reshape(dataBOLD.catmvoltc, nVox, nt)', S(iC).rgrsMION_resample','rows','complete');

    S(iC).mapR = reshape(Rvals, [nx, ny, nz]).*(-1); % because of MION
    S(iC).mapP = reshape(Pvals, [nx, ny, nz]);
    fprintf(1, '%d/%d: cell %s\n', iC, length(S), S(iC).cellID)
        toc
end



%% Correlation across trials
% sliding window correlation across trials

sizeW = 5000; % window size in msec
fprintf(1, 'Compute correlation across trials:: \n')
for iC = 1:length(S)
    tic
    fprintf(1, '%d/%d: cell %s \n', iC, length(S), S(iC).cellID)
    [C, Cparams] = computeSlidingWindowCorrAcrossTrials(S(iC).matsdf, sizeW);
    S(iC).corrTrials = C;
    S(iC).corrTrialsParams = Cparams;
    toc
end

%
stepW = S(iC).corrTrialsParams.stepW;
sizeW = S(iC).corrTrialsParams.sizeW;
numW = S(iC).corrTrialsParams.numW;
axisWindow = sizeW/2:stepW:stepW*(numW-1)+sizeW/2;
% iM=1;

nRowSubplot = length(setMovIDs)/3+1;
figure;
for iC=1:length(S)
    fprintf(1, 'Channel %s (%d/%d)', S(iC).cellID, iC, length(S));
    sp1=subplot(nRowSubplot,3,[1 2 3]);
    plot(S(iC).catmnsdf);
    title(sp1, sprintf('Channel %s: SDF averaged across trials: concatenated for %d movies', S(iC).cellID, length(setMovIDs)))
    axis tight;
    for iM=1:length(S(iC).movID)
    sp2=subplot(nRowSubplot,3,iM+3);    
    plot(squeeze(S(iC).matsdf{iM}(:,:))) %plot(squeeze(S_whole(iC).matsdf{iM}(:,:))) %plot(squeeze(S(iC).matsdf(:,:,iM)))
    title(sp2, sprintf('Movie #%d', S(iC).movID(iM)))
    axis tight
    end
%     sp3=subplot(5,3,[13 14 15]);    
%     plot(axisWindow, mean(S(iC).corrTrials(iM).R), 'o')
%     title(sp3,'Sliding-window correlation across trials')
    input('')
end

% save SDF for individual trials separately (it's too large)
tempS = S;
Names = fieldnames(S);
tempS = rmfield(tempS, Names([1:2, 4:11], 1));

save(fullfile(dirDataNeural, saveFileName_matSDF), 'tempS','-v7.3')

% save regressors & BOLD 
S = rmfield(S, 'matsdf');
save(fullfile(dirData, saveFileName), 'S', 'dataBOLD', '-v7.3')



% speed4=real(log(speed_rgr));          % compress with a log function
% speed4(isinf(speed4))=0;              % make sure we have no INF numbers
% speed5=doConv(speed4,MION_k)';        % apply MION function to raw data
% 
% dmot5=doConv(dmot4,MION_k)';        % apply MION function to raw data
% dmot6=doConv(dmot5,smooth_k)';      % apply a 2 standard deviation smoothing kernel

  