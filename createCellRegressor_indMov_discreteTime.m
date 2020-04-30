function FR_dT = createCellRegressor_indMov_discreteTime(dirDataNeural, setCellIDs, setMovIDs, sizeTimeBin_sec) % , indDataMov) %, flagConcat)

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
d_n = dir(fullfile(dirDataNeural, '*sig*.mat'));
d_n = d_n(cellfun(@length, {d_n(:).name})<30);
[listMatSUFile{1:length(d_n), 1}] = deal(d_n.name);

% list of channel IDs and movie IDs of each file
listSU_all = regexp(cat(2,d_n.name), '(?<=sig)\w*', 'match')'; % list of channels
listMov_all = regexp(cat(2,d_n.name), '\d*(?=sig)', 'match')'; % list of movies

% Since not every channel has data for every movie, 
% the "data presence matrix (Channels x MovieIDs)" is generated here
% for which channel has which movie's data: 1 for data, 0 for no data
listSUchannelID = unique(listSU_all);

indDataMov=[];
for iCh=1:length(listSUchannelID) % go throucgh channel-by-channel     
     tempListMov{iCh} = regexp([listMatSUFile{strcmp(listSUchannelID{iCh}, listSU_all)}], '\d*(?=sig)', 'match');
     % Get indices for common movies across cells
     indDataMov(iCh, :) = ismember(setMovIDs, str2num(char(tempListMov{iCh}))'); % data presence matrix (Channels x MovieIDs) 
end




% %cellIDs = {'003a';'007a';'013a';'013b';'014a'};
% cellIDs = {'003a'};
% %movIDs  = [1 2 3 7 8 9];
% movIDs  = [7 8 9];

% % Convolve MION function with mean SDF, then resample
% fMRI_TR_sec = 2.4;
% timeResNeural_sec = 0.001; %f


FR_dT = [];

for iChan = 1:length(setCellIDs)
  
  cur_cellID = setCellIDs{iChan};
  indCurChan = find(strcmp(cur_cellID, listSUchannelID)>0);
%   count = 0;
  for iMov = 1:length(setMovIDs) % should concatenate all those movies in order
      
      cur_movID = setMovIDs(iMov);
      
      if ~indDataMov(indCurChan, iMov), 
          fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(indCurChan, iMov))
          fprintf(1, 'Skip to the next movie/cell \n')
          continue; 
      end
      
      % Find a relevant file
      filename = char(listMatSUFile((strcmp(cur_cellID, listSU_all)+strcmp(num2str(cur_movID), listMov_all))==2));
      fprintf(1, 'cell: %s, movie: %d, data: %d \n', cur_cellID, cur_movID, indDataMov(indCurChan, iMov))
      fprintf(1, 'filename: %s\n', filename)

      FR_dT(iChan, iMov).cellID = cur_cellID;
      FR_dT(iChan, iMov).movID = cur_movID;
      FR_dT(iChan, iMov).SUfilename = filename;
      [FR_dT(iChan, iMov).matFR{1}, FR_dT(iChan, iMov).mnFR] = computeMeanFR(dirDataNeural, filename, sizeTimeBin_sec); % averaged SDF across trials

%       [S(iChan, iMov).matsdf{1}, S(iChan, iMov).mnsdf] = computeMeanSDF(dirDataNeural, filename); % averaged SDF across trials
%       [S(iChan, iMov).rgrsMION, S(iChan, iMov).rgrsMION_resample] = computeMIONrgrs(SDF_dT(iChan, iMov).mnsdf, timeResNeural_sec, fMRI_TR_sec);

      
  end
      
  
end
    
function [matFR, mnFR]= computeMeanFR(dirDataNeural, filename, sizeTimeBin_sec)
  
% global dat flagLocal
% 
% if flagLocal
%     addpath('/Volumes/share/UCNI_Library/matlab_utils/');
%     dirData = '/Volumes/USRlab/data/parksh/'; %/rmov/';
% else
%     addpath('/einstein0/share/UCNI_Library/matlab_utils/');
%     dirData = '/einstein0/USRlab/data/parksh/'; %/rmov/';
% end
load(fullfile(dirDataNeural, filename))

nTrial = length(dat.s);

% prop = 'Color';
% val = 'k';

switch lower(dat.h.units)
    case 'sec'
        win = [0 300]; % seconds or milliseconds
    case 'ms'
        win = [0 300*1000];
end

spikes = {};
for t=1:nTrial %[1 2 5 6 7 8]
    ts = dat.s{t};
    spikes{t} = ts((ts>=win(1)) & ts<=win(2));
end

% % sdf sample period=;/fghipr
% sdf_samp = 100;
% 

matFR = [];
for i=1:length(spikes)
  matFR(:,i) = getFR(spikes{i}, win, sizeTimeBin_sec);
end

mnFR = mean(matFR,2);


function fr_out = getFR(ts, win, sizeTimeBin_sec)
upsamp = 1000;
%   TR = 2.4;
%   MION_k = gampdf(-3*TR*upsamp:3*TR*upsamp,TR,2*TR);
%   MION_k = MION_k./sum(MION_k);

if max(win)<1000 % if window is in seconds, we need to upsample it
    % it's not the exact way to do this, but just for now..
    timeline = zeros(upsamp*win(2)-upsamp*win(1),1);
    its = (ceil(ts*upsamp));
%     timeline(its) = 1;
    
else % units are milliseconds
    timeline = zeros(win(2)-win(1),1);
    its = ceil(ts);
%     timeline(its) = 1;
    
end

sizeTimeBin_ms = sizeTimeBin_sec*1000;
nBin = round(length(timeline)/sizeTimeBin_ms);

timeBin = [0:nBin].*sizeTimeBin_ms;
% if timeBin(end) < length(timeline) % if there is a remaining window 
%     timeBin = cat(2, timeBin, length(timeline));
% end

fr_out = zeros(nBin,1);
for iBin = 1:nBin
    fr_out(iBin,1) = sum(its >= timeBin(iBin)+1 & its < timeBin(iBin+1));
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

  