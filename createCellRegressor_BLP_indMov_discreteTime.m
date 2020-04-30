function [BLPRGR] = createCellRegressor_BLP_indMov_discreteTime(dirDataNeural, sizeTimeBin_sec, flagSave) %setCellIDs, setMovIDs,  %indDataMov) %, flagConcat)
%
% Create fMRI regressor with Band Limited Power (BLP) of LFP signals
% For each movie, each band separately saved in structure format
%
% 2015/02/05 SHP
%
% 2015/02/09 SHP
% Replace data points exceeding 5*std amplitude with NaNs
%
% 2015/03/19 SHP
% 1. Sort the movie ID from small to big, so that the resulting struct is
%   ordered
% 2. Add the case when the time unit in file is 'sec' instead of 'ms'
%
% 2015/03/25 SHP
% 1. flagSave: 1 for save the result (default: no save)

if nargin < 3
    flagSave = 0;
end

% clear global S
% global S

% global flagLocal
%
% if flagLocal
%     addpath('/Volumes/share/UCNI_Library/matlab_utils/');
%     dirData = '/Volumes/USRlab/data/parksh/'; %/rmov/';
% else
%     addpath('/einstein0/share/UCNI_Library/matlab_utils/');
%     dirData = '/einstein0/USRlab/data/parksh/'; %/rmov/';
% end

% dirData = '/einstein0/USRlab/data/parksh/'; %/rmov/';
% W = what(dirDataNeural);
d_n = dir(fullfile(dirDataNeural, '*LFP.mat')); %'*sig*.mat'));
[listLFPFile{1:length(d_n), 1}] = deal(d_n.name);

% list of channel IDs and movie IDs of each file
listMov_all = regexp(cat(2,d_n.name), '\d*(?=LFP)', 'match')'; % list of movies
% listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')'; % list of channels


% % Since not every channel has data for every movie,
% % the "data presence matrix (Channels x MovieIDs)" is generated here
% % for which channel has which movie's data: 1 for data, 0 for no data
% listSUchannelID = unique(listSU_all);
% indDataMov=[];
% for iCh=1:length(listSUchannelID) % go throucgh channel-by-channel
%      tempListMov{iCh} = regexp([listLFPFile{strcmp(listSUchannelID{iCh}, listSU_all)}], '\d*(?=sig)', 'match');
%      % Get indices for common movies across cells
%      indDataMov(iCh, :) = ismember(setMovIDs, str2num(char(tempListMov{iCh}))'); % data presence matrix (Channels x MovieIDs)
% end

% %cellIDs = {'003a';'007a';'013a';'013b';'014a'};
% cellIDs = {'003a'};
% %movIDs  = [1 2 3 7 8 9];
% movIDs  = [7 8 9];

% Convolve MION function with mean SDF, then resample
fMRI_TR_sec = 2.4;
timeResNeural_sec = 0.005; %5ms  %0.001; %f


setMovIDs = sort(str2double(listMov_all)); % for all movies
BLPRGR = struct();
% for iChan = 1:length(setCellIDs)
%
%   cur_cellID = setCellIDs{iChan};
%   count = 0;
for iMov = 1:length(setMovIDs) % should concatenate all those movies in order
    
    cur_movID = setMovIDs(iMov);
    
    %       if ~indDataMov(iChan, iMov),
    %           fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(iChan, iMov))
    %           fprintf(1, 'Skip to the next movie/cell \n')
    %           continue;
    %       end
    
    % Find a relevant file
    filename = char(listLFPFile(strcmp(num2str(cur_movID), listMov_all)));
    fprintf(1, 'movie: %d, filename: %s\n', cur_movID, filename)
    
    
    %       BLP_dT(iChan, iMov).cellID = cur_cellID;
    BLPRGR(iMov).movID = cur_movID;
    BLPRGR(iMov).filename = filename;
    [BLPRGR(iMov).matBLP, BLPRGR(iMov).meanBLP, BLPRGR(iMov).blpfreq] = computeMeanFR(dirDataNeural, filename, sizeTimeBin_sec); % averaged SDF across trials
    
    %       BLP_dt(iMov).blpfreq = dat.h.blpfreq;
    
    %       [S(iChan, iMov).matsdf{1}, S(iChan, iMov).mnsdf] = computeMeanSDF(dirDataNeural, filename); % averaged SDF across trials
    %       [S(iChan, iMov).rgrsMION, S(iChan, iMov).rgrsMION_resample] = computeMIONrgrs(SDF_dT(iChan, iMov).mnsdf, timeResNeural_sec, fMRI_TR_sec);
    
end

if flagSave
    [dirDataHome, nameSubj, c] = fileparts(dirDataNeural);
    saveFileName = sprintf('%s_movieTS_BLPLFP_indMov.mat', nameSubj);
    save(fullfile(dirDataNeural, saveFileName), 'BLPRGR')
    fprintf(1, '\n Data is saved as "%s" \n', saveFileName)
end


% end

function [matBLPFR, meanBLP, blpfreq]= computeMeanFR(dirDataNeural, filename, sizeTimeBin_sec)


load(fullfile(dirDataNeural, filename))

nTrial = size(dat.blp,1);
nBand = size(dat.blp,3);

% prop = 'Color';
% val = 'k';

switch lower(dat.h.units)
    case 'sec'
        win = [0 300]; % seconds or milliseconds
        time = dat.t(dat.t>=win(1) & dat.t<=win(2));
        time = time*1000;
    case 'ms'
        win = [0 300*1000]; % individual 5-min movie in milliseconds
        time = dat.t(dat.t>=win(1) & dat.t<=win(2));
end
% tloc = find(dat.t>=win(1) & dat.t<=win(2));

meanBLP = {};
critStd = 5; % criterion for extreme value (in std)
for iBP = 1:nBand % for each frequency band
    %     blp = {};
    
    % For each trial, change the extreme values to NaNs
    matBLP = squeeze(dat.blp(:,:, iBP))'; % time x trials
    matBLP_zscore = zscore(matBLP, 0, 1); % zsocred for each trial
    
    validLoc = abs(matBLP_zscore)<critStd;
    
    matBLP_ext = NaN(size(matBLP));
    matBLP_ext(validLoc) = matBLP(validLoc); % original value in microvolts
    
    
    %     for t=1:nTrial %[1 2 5 6 7 8]
    %         ts = squeeze(dat.blp(t, :, iBP));
    %         blp{t} = ts(tloc);
    % %         ts = dat.s{t};
    % %         blp{t} = ts((ts>=win(1)) & ts<=win(2));
    %     end
    
    % % sdf sample period=;/fghipr
    % sdf_samp = 100;
    %
    
    tmpmatBLP = [];
    for iTR=1:nTrial
        tmpmatBLP(:,iTR) = getFR(matBLP_ext(:,iTR), win, time, sizeTimeBin_sec); %getFR(blp{iTR}, win, time, sizeTimeBin_sec);
    end
    
    mnFR = nanmean(tmpmatBLP, 2); %mean(matFR,2);
    meanBLP{iBP} = mnFR;
    matBLPFR{iBP} = tmpmatBLP;
    
end

blpfreq = dat.h.blpfreq;



function fr_out = getFR(ts, win, time, sizeTimeBin_sec)
upsamp = 1000;
%   TR = 2.4;
%   MION_k = gampdf(-3*TR*upsamp:3*TR*upsamp,TR,2*TR);
%   MION_k = MION_k./sum(MION_k);



if max(win)<1000 % if window is in seconds, we need to upsample it
    % it's not the exact way to do this, but just for now..
    timeline = zeros(upsamp*win(2)-upsamp*win(1),1);
    %     its = (ceil(ts*upsamp));
    %     timeline(its) = 1;
    
else % units are milliseconds
    timeline = zeros(win(2)-win(1),1);
    %     its = ceil(ts);
    %     timeline(its) = 1;
    
end

sizeTimeBin_ms = sizeTimeBin_sec*1000;
nBin = length(timeline)/sizeTimeBin_ms;

timeBin = [0:nBin].*sizeTimeBin_ms;

fr_out = zeros(nBin,1);
for iBin = 1:nBin
    fr_out(iBin,1) = nanmean(ts(time >= timeBin(iBin)+1 & time < timeBin(iBin+1))); %sum(its >= timeBin(iBin)+1 & its < timeBin(iBin+1));
end




%% Compute Spike Density Function and MION convolution
% function [matSdf, mnsdf]= computeMeanSDF(dirDataNeural, filename)
%
% % global dat flagLocal
% %
% % if flagLocal
% %     addpath('/Volumes/share/UCNI_Library/matlab_utils/');
% %     dirData = '/Volumes/USRlab/data/parksh/'; %/rmov/';
% % else
% %     addpath('/einstein0/share/UCNI_Library/matlab_utils/');
% %     dirData = '/einstein0/USRlab/data/parksh/'; %/rmov/';
% % end
% load(fullfile(dirDataNeural, filename))
%
% nTrial = length(dat.s);
%
% % prop = 'Color';
% % val = 'k';
%
% switch lower(dat.h.units)
%     case 'sec'
%         win = [0 300]; % seconds or milliseconds
%     case 'ms'
%         win = [0 300*1000];
% end
%
% spikes = {};
% for t=1:nTrial %[1 2 5 6 7 8]
%     ts = dat.s{t};
%     spikes{t} = ts((ts>=win(1)) & ts<=win(2));
% end
%
% % % sdf sample period=;/fghipr
% % sdf_samp = 100;
% %
%
% matSdf = [];
% for i=1:length(spikes)
%   matSdf(:,i) = getSDF(spikes{i},win);
% end
%
% mnsdf = mean(matSdf,2);
%
% function sdf_out = getSDF(ts, win)
%
% sd = 1000; %2000; %1000; %2000;
% k  = normpdf(-3*sd:+3*sd,0,sd);
% k = k./sum(k);
% upsamp = 1000;
% %   TR = 2.4;
% %   MION_k = gampdf(-3*TR*upsamp:3*TR*upsamp,TR,2*TR);
% %   MION_k = MION_k./sum(MION_k);
%
% if max(win)<1000 % if window is in seconds, we need to upsample it
%     % it's not the exact way to do this, but just for now..
%     timeline = zeros(upsamp*win(2)-upsamp*win(1),1);
%     its = (ceil(ts*upsamp));
%     timeline(its) = 1;
%
% else % units are milliseconds
%     timeline = zeros(win(2)-win(1),1);
%     its = ceil(ts);
%     timeline(its) = 1;
%
% end
%
% % convolve, scale to spikes/sec, and decimate
% sdf_out = doConv(timeline,k).*upsamp; %doConv(timeline, MION_k).*upsamp; %


%
% function [rgrsMION, rgrsMION_resample] = computeMIONrgrs(catmnsdf, timeResNeural_sec, fMRI_TR_sec)
% % % set up MION function for convolution
% % tvals = [(-20*dataBOLD.TR):dataBOLD.TR:(20*dataBOLD.TR)];
% % % This is just made-up
% % MION_k = gampdf(tvals,dataBOLD.TR,2*dataBOLD.TR);
% % MION_k = MION_k./sum(MION_k);
%
%
% % another MION function
% lengthMovie_sec = 300; %length(catmnsdf)/1000/9; %3;
% t = 0:fMRI_TR_sec:lengthMovie_sec;
% h = 16.4486 * ( -0.184/ 1.5 * exp(-t/ 1.5)...
% +0.330/ 4.5 * exp(-t/ 4.5)...
% +0.670/13.5 * exp(-t/13.5) );
%
% % figure;
% % plot(t, h, 'o-')
%
% % 'MION(d)'
% % = 1 parameter block stimulus of duration 'd',
% % intended to model the response of MION.
% % The zero-duration impulse response 'MION(0)' is
% % h(t) = 16.4486 * ( -0.184/ 1.5 * exp(-t/ 1.5)
% % +0.330/ 4.5 * exp(-t/ 4.5)
% % +0.670/13.5 * exp(-t/13.5) )
% % which is adapted from the paper
% % FP Leite, et al.  NeuroImage 16:283-294 (2002)
% % http://dx.doi.org/10.1006/nimg.2002.1110
% % ** Note that this is a positive function, but MION
% % produces a negative response to activation, so the
% % beta and t-statistic for MION are usually negative.
% % ***** If you want a negative MION function (so you get
% % a positive beta), use the name 'MIONN' instead.
% % ** After convolution with a square wave 'd' seconds
% % long, the resulting single-trial waveform is
% % scaled to have magnitude 1.  For example, try
% %     this fun command to compare BLOCK and MION:
% %     3dDeconvolve -nodata 300 1 -polort -1 -num_stimts 2   \
% %     -stim_times 1 '1D: 10 150' 'MION(70)'    \
% %     -stim_times 2 '1D: 10 150' 'BLOCK(70,1)' \
% %     -x1D stdout: | 1dplot -stdin -one -thick
% %     You will see that the MION curve rises and falls
% %     much more slowly than the BLOCK curve.
% %     ==> ** Note that 'MION(d)' is already convolved with a
% %     square wave of duration 'd' seconds.  Do not
% %     convolve it again by putting in multiple closely
% %     spaced stimulus times (this mistake has been made)!
% %     ** Scaling the single-trial waveform to have magnitude
% %     1 means that trials with different durations 'd'
% %     will have the same magnitude for their regression
% %     models.
%
%
% % % set up standard smoothing kernel
% % smooth_k = normpdf([-10:10],0,2);
% % smooth_k = smooth_k/sum(smooth_k);
%
% rgrsMION = doConv(catmnsdf, h); %MION_k);
% % rgrsMION_sm = doConv(rgrsMION, smooth_k);
%
% rgrsMION_resample = resample(rgrsMION, timeResNeural_sec*1000, fMRI_TR_sec*1000); % inputs need to be integer
%

