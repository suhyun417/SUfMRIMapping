function [S, dataBOLD] = getCorrMap_indMov(nameSubjNeural, nameSubjBOLD,saveFileName_matSDF, saveFileName) %(nameSubjNeural, nameSubjBOLD, saveFileName_matSDF, saveFileName)

% % Make individual
% nameSubjNeural = input('Subject name for single unit data? (e.g. Tor, Rho)', 's');
% nameSubjBOLD = input(

%% Environmental setup
% global flagLocal

flagLocal = 0; %2; %0; %1; % 0 for server (Pauli), 1 for local machine, 2 for home

switch flagLocal
    case 0 % Pauli
    addpath('/einstein0/share/UCNI_Library/matlab_utils/');
    dirData = '/einstein0/USRlab/data/parksh/'; %/rmov/';
    case 1 % SP's desktop
    addpath('/Volumes/share/UCNI_Library/matlab_utils/');
    dirData = '/Volumes/USRlab/data/parksh/'; %/rmov/';
%     case 2 % home
%     addpath('/Users/shyuny/Documents/MATLAB/expToolbox/BlockAna/'); %'/einstein0/share/UCNI_Library/matlab_utils/');
%     addpath('/Users/shyuny/Documents/MATLAB/expToolbox/BlockAna/matlab_utils');
%     dirData = '/Volumes/HONGKONG/USRlab/data/parksh/'; %'/Users/shyuny/ResearchProjects/0NeuralBOLDCorr/data/';
end
% W = what(dirData);

% Subject name and directory
% nameSubjNeural = 'Tor'; %'Rho';
% nameSubjBOLD = 'Art';
dirDataNeural = [dirData, nameSubjNeural, '/'];
dirDataBOLD = [dirData, nameSubjBOLD, '/'];


%%
% Get the movie list from single unit data
d_n = dir(fullfile(dirDataNeural, '*sig*.mat'));
[listMatSUFile{1:length(d_n), 1}] = deal(d_n.name); % Get the list of single unit files 
% listMatSUFile = regexp(cat(2, d_n.name), '\w*.mat', 'match')'; %regexp(cat(2, d_n.name), '(?<=mat)\w*', 'match')';
% listSUchannelID = unique(regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match'));
listMovSU = unique(regexp(cat(2,d_n.name), '\d*(?=sig)', 'match')); % Movie IDs of the single units in cellstr

% Get the movie list from the fMRI data
fprintf(1, '\nPreparing BOLD time series...\n');

d_b = dir(fullfile(dirDataBOLD, '*BOLD.mat')); %W = what(dirDataBOLD);
load(fullfile(dirDataBOLD, d_b(1).name)); % for now it's only one file
listMovBOLD = MCD.unimov; % Movie IDs of the current BOLD



% Select common movie IDs between single units and BOLD (or define the movies)
setMovIDs = input('Movie IDs (e.g. [1 2 3])? (default: common movies in single units and fMRI)');
if isempty(setMovIDs)
    setMovIDs = intersect(listMovBOLD, str2num(char(listMovSU))');
end
% setMovIDs = [13 14 15]; %intersect(listMovBOLD, str2num(char(listMovSU))');
[~, locMovBOLD] = intersect(listMovBOLD, setMovIDs);
% [setMovIDs, locMovBOLD] = intersect(listMovBOLD, str2num(char(listMovSU))'); 

% flagConcat = input('Do you want to concatenate timeseries across movies? (1:yes, 0:no):');

% Get the list of single units & their movie
listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')';
listSUchannelID = unique(listSU_all);

% indDataMov=[];
% for iCh=1:length(listSUchannelID) % go throucgh channel-by-channel     
%      tempListMov{iCh} = regexp([listMatSUFile{strcmp(listSUchannelID{iCh}, listSU_all)}], '\d*(?=sig)', 'match');
%      % Get indices for common movies across cells
%      indDataMov(iCh, :) = ismember(setMovIDs, str2num(char(tempListMov{iCh}))'); % data presence matrix (Channels x MovieIDs) 
% end

% % For now, let's do it with cells that have data for every movies
% setCellIDs = listSUchannelID(sum(indMov,2)==length(setMovIDs));
setCellIDs = listSUchannelID;

% setCellIDs = {'006a'; '009a'; '010a'; '012a'; '013a'}; % Rhombus %{'003a';'005a';'007a';'009a';'013b';'014a'}; % Sig: '013a' has very low spikes
% setMovIDs  = [7 8 9 10 11 12]; % Rhom %[7 8 9]; % Sig % [1 2 3 7 8 9];





%% First, split fMRI responses into each individual movie

% fprintf(1, 'Number of movies in original time series: %d\n', length(setMovIDs))

dataBOLD = MCD; % Get necessary params for further DSP in BlockAna: for now, just all the params
dataBOLD = rmfield(dataBOLD, 'catmvoltc'); % remove the previous time series

lengthIndMovie_TR = dataBOLD.max_tr;
% lengthIndMovie_sec = lengthIndMovie_TR*dataBOLD.TR; % 300;

% Get the necessary part of time courses, corresponding to our (common) movie set
% if ~issorted(setMovIDs) % just to be sure
%     setMovIDs = sort(setMovIDs);
% end

% countTR = 0; tc=[];
for iMov = 1:length(locMovBOLD) % concatenate across different movies
    fprintf(1, 'Movie ID: %d (%d/%d), Location in original BOLD: %d\n', setMovIDs(iMov), iMov, length(setMovIDs), locMovBOLD(iMov));
%     tc(:,:,:,countTR+1:countTR+lengthIndMovie_TR) = MCD.catmvoltc(:,:,:,lengthIndMovie_TR*(locMovBOLD(iMov)-1)+1:lengthIndMovie_TR*locMovBOLD(iMov)); % should be optimized for different set of movies
    dataBOLD.catmvoltc{iMov} = MCD.catmvoltc(:,:,:,lengthIndMovie_TR*(locMovBOLD(iMov)-1)+1:lengthIndMovie_TR*locMovBOLD(iMov)); % should be optimized for different set of movies
    %     countTR = size(tc, 4);
end
% dataBOLD.catmvoltc = tc;
% fprintf(1, 'Time course for %d movies has been saved in dataBOLD.tc\n', length(setMovIDs));
% fprintf(1, 'The number of TRs of dataBOLD.catmvoltc should be %d and is %d\n', ...
%     lengthIndMovie_TR*length(setMovIDs), size(dataBOLD.catmvoltc, 4));

dataBOLD.unimov = setMovIDs;

clear MCD %tc %SigData







%% Get mean SDF and regressors for each single unit
S = createCellRegressor_indMov(dirDataNeural, setCellIDs, setMovIDs) %, indDataMov); %, flagConcat); 
% -- S(iChan, iMov): for each cell & each movie (S.cellID), 
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
for iC = 1:size(S,1) %length(S)
  
    for iM = 1:size(S,2)
        fprintf(1, 'Channel %s (%d/%d), Movie %d', S(iC).cellID, iC, size(S,1), setMovIDs(iM));
        figure(fig_mnsdf);
        %   sp(iC) = subplot(length(S),1,iC);
        plot(S(iC, iM).mnsdf); %plot(S(iC).catmnsdf);
        axis tight;
        %   title(sprintf('cell: %s, movies: [%d %d %d]', S(iC).cellID, setMovIDs))
        
        figure(fig_rgrsMIONresample);
        %   sp2(iC) = subplot(length(S),1,iC);
        plot(S(iC, iM).rgrsMION_resample); %plot(S(iC).rgrsMION_resample);
        axis tight;
        %   title(sprintf('cell: %s, movies: [%d %d %d]', S(iC).cellID, setMovIDs))
        input('')
    end
    
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
% figure;
% set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto', 'Name', 'MION regressors and BOLD response from one example voxel')
% sp3=[];
for iC = 1:size(S,1) %length(S)
%     fprintf(1, 'Channel %s (%d/%d)', S(iC, iM).cellID, iC, length(S));
    for iM = 1:size(S,2) %length(setMovIDs)         
        if ~isempty(S(iC, iM).rgrsMION_resample)
            fprintf(1, 'Channel %s (%d/%d), Movie %d \n', S(iC).cellID, iC, size(S,1), setMovIDs(iM));
            S(iC, iM).rgrsMION_resample(1:onsetResp_numTR) = NaN; %S(iC, iM).rgrsMION_resample(lengthIndMovie_TR*(iM-1)+1:lengthIndMovie_TR*(iM-1)+onsetResp_numTR) = NaN;
%             plot(S(iC, iM).rgrsMION_resample);
%             input('')
        end
    end
%     sp3(iC) = subplot(length(S)+1,1,iC);
%     plot(S(iC).rgrsMION_resample);     
%     title(sprintf('cell: %s, movies: [%d %d %d]', S(iC).cellID, setMovIDs))
%     input('')
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
[nx, ny, nz, nt] = size(dataBOLD.mvoltc{1});
nVox = nx*ny*nz;

% % 1. gamma pdf 
TR=2.4;
k = gampdf([-40:TR:40],4,2);
% % 2. kernel from AFNI
% taxis = 0:2.4:50;
% k = 16.4486 * ( -0.184/ 1.5 * exp(-taxis/ 1.5)...
% +0.330/ 4.5 * exp(-taxis/ 4.5)...
% +0.670/13.5 * exp(-taxis/13.5) );
% k = cat(2, zeros(1, length(k )-1), k );

fprintf(1, 'Compute correlation using regressors:: \n')
for iC = 1:size(S,1) %length(S)  
    for iM = 1:size(S,2)
        
        if ~isempty(S(iC, iM).mnFR)
            tic
            
            neuralrgrs = S(iC, iM).mnFR;
            neuralrgrs = neuralrgrs-mean(neuralrgrs);
            neuralrgrs = doConv(neuralrgrs,k);%conv(neuralrgrs,k,'same');
            
            [Rvals, Pvals] = corr(reshape(dataBOLD.mvoltc{iM}, nVox, nt)', neuralrgrs','rows','complete');
            
            matCorr(iC, iM).mapR = reshape(Rvals, [nx, ny, nz]).*(-1); % because of MION
            matCorr(iC, iM).mapP = reshape(Pvals, [nx, ny, nz]);
            fprintf(1, 'cell %s (%d/%d), movie: %d \n', S(iC, iM).cellID, iC, size(S,1), iM) %fprintf(1, '%d/%d: cell %s\n', iC, length(S), S(iC, iM).cellID)
            toc
        end
    end
end



%% Correlation across trials
% sliding window correlation across trials

sizeW = 5000; % window size in msec
fprintf(1, 'Compute correlation across trials:: \n')
for iC = 1:size(S,1) %length(S)
    for iM = 1:size(S,2)
        if ~isempty(S(iC, iM).rgrsMION_resample)
            tic
            fprintf(1, 'cell %s (%d/%d), movie: %d \n', S(iC, iM).cellID, iC, length(S), iM) %fprintf(1, '%d/%d: cell %s \n', iC, length(S), S(iC, iM).cellID)
            [C, Cparams] = computeSlidingWindowCorrAcrossTrials(S(iC, iM).matsdf, sizeW);
            S(iC, iM).corrTrials = C;
            S(iC, iM).corrTrialsParams = Cparams;
            toc
        end
    end
end

% % Plotting the results
% stepW = S(iC, iM).corrTrialsParams.stepW;
% sizeW = S(iC, iM).corrTrialsParams.sizeW;
% numW = S(iC, iM).corrTrialsParams.numW;
% axisWindow = sizeW/2:stepW:stepW*(numW-1)+sizeW/2;
% % iM=1;
% 
% nRowSubplot = length(setMovIDs)/3+1;
% figure;
% for iC=1:length(S)
%     fprintf(1, 'Channel %s (%d/%d)', S(iC, iM).cellID, iC, length(S));
%     sp1=subplot(nRowSubplot,3,[1 2 3]);
%     plot(S(iC, iM).catmnsdf);
%     title(sp1, sprintf('Channel %s: SDF averaged across trials: concatenated for %d movies', S(iC, iM).cellID, length(setMovIDs)))
%     axis tight;
%     for iM=1:length(S(iC, iM).movID)
%     sp2=subplot(nRowSubplot,3,iM+3);    
%     plot(squeeze(S(iC, iM).matsdf{iM}(:,:))) %plot(squeeze(S_whole(iC).matsdf{iM}(:,:))) %plot(squeeze(S(iC, iM).matsdf(:,:,iM)))
%     title(sp2, sprintf('Movie #%d', S(iC, iM).movID(iM)))
%     axis tight
%     end
% %     sp3=subplot(5,3,[13 14 15]);    
% %     plot(axisWindow, mean(S(iC).corrTrials(iM).R), 'o')
% %     title(sp3,'Sliding-window correlation across trials')
%     input('')
% end

% save SDF for individual trials separately (it's too large)
tempS = S;
Names = fieldnames(S);
tempS = rmfield(tempS, Names([1:3, 5:9], 1)); %rmfield(tempS, 'matsdf'); %Names([1:2, 4:11], 1));

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

  