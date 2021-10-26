function computeCorr
% This function converts LFP input into band-limited power bands and then 
% calculates the correlations between all channels for these power bands.
% USAGE:    computerCorr
% INPUT:    LFP is the 2-dimensional LFP data, such that the first 
%           dimension represents channel and the second dimension represents 
%           time. fs is the sampling rate; you should have your data already 
%           downsampled.
%           resFile is the name of the resulting BLP structure
% OUTPUT:   BLP (structure with 3-dim matrices)
%           BLP.MA.blp_dat w/ [channel|time|bandwidth] of Mystery Area
%           BLP.V1.blp_dat w/ [channel|time|bandwidth] of V1
%           BLP.fs_BLP
%           BLP.frange

date = '21-05-08';  date2 = '210508';  sess = '_0001';
fs = 1000;

cd(['/einstein0/USRlab/projects/scholvinckm/data/Varia/outside scanner/' date '/Matlab']);
load([date2 sess '_LFP_1000']);

LFP_MA = LFP_MA';
LFP_V1 = LFP_V1';

nyq = fs/2;
numchan = size(LFP_MA,2);
frange = getFrequencyRanges;
nbw    = size(frange,2);

% for MA
for chan = 1:numchan

   rdat = squeeze(LFP_MA(:,chan));
   filt_dat = zeros(size(rdat,1),nbw);
   blp_dat = zeros(size(rdat,1),nbw);
       
   for bw = 1:nbw
        fr                    = frange(:,bw);
        hpc                   = fr(1);                      % Hz
        lpc                   = fr(2);                      % Hz
        hWn                   = hpc/nyq;                    % between 0.0 and 1.0
        lWn                   = lpc/nyq;                    % between 0.0 and 1.0
        [lb,la]               = cheby1(2,0.8,[hWn lWn]);    % chebyshev bandpass
        fdat                  = filtfilt(lb,la,rdat);       % bandpass
        filt_dat(:,bw)        = fdat;
        blp_dat(:,bw)         = abs(fdat);                  % "power"
        meancorr_blp(:,bw)    = abs(fdat) - mean(abs(fdat));% mean correct blp
   end

    BLP.MA.blp_dat(chan,:,:)   = meancorr_blp;              % channel x timepoint x bandwidth
    BLP.fs_BLP   = fs;
    BLP.frange   = frange;
    
end

% for V1
for chan = 1:numchan

   rdat = squeeze(LFP_V1(:,chan));
   filt_dat = zeros(size(rdat,1),nbw);
   blp_dat = zeros(size(rdat,1),nbw);
       
   for bw = 1:nbw
        fr                    = frange(:,bw);
        hpc                   = fr(1);                      % Hz
        lpc                   = fr(2);                      % Hz
        hWn                   = hpc/nyq;                    % between 0.0 and 1.0
        lWn                   = lpc/nyq;                    % between 0.0 and 1.0
        [lb,la]               = cheby1(2,0.8,[hWn lWn]);    % chebyshev bandpass
        fdat                  = filtfilt(lb,la,rdat);       % bandpass
        filt_dat(:,bw)        = fdat;
        blp_dat(:,bw)         = abs(fdat);                  % "power"
        meancorr_blp(:,bw)    = abs(fdat) - mean(abs(fdat));% mean correct blp
   end

    BLP.V1.blp_dat(chan,:,:)   = meancorr_blp; % channel x timepoint x bandwidth
    
end

BLP_MA = permute(BLP.MA.blp_dat(:,:,:),[2 1 3]); % need permute for corrcoef to work
BLP_V1 = permute(BLP.V1.blp_dat(:,:,:),[2 1 3]);

for bw = 1:nbw
    
    % compute correlations
    corrMap = zeros(2*numchan,2*numchan);
    % MA with MA
    corrMap(1:numchan,1:numchan) = corrcoef(BLP_MA(:,1:numchan,bw));
    % MA with V1
    for j=1:numchan
        for i=1:numchan
            corr = corrcoef(BLP_MA(i,:,bw),BLP_V1(j,:,bw));
            corrMap(i,j+numchan) = corr(1,2);
        end
    end
    % V1 with MA
    for j = 1:numchan
        for i=1:numchan    
            corr = corrcoef(BLP_V1(i,:,bw),BLP_MA(j,:,bw));
            corrMap(i+numchan,j) = corr(1,2);
        end
    end
    % V1 with V1
    corrMap(numchan+1:2*numchan,numchan+1:2*numchan) = corrcoef(BLP_V1(:,1:numchan,bw));
    
    % plot correlations
    figure(bw); imagesc(corrMap); colormap(jet(256)); colorbar;
    
    Corr(bw).corrmap = corrMap;
    
end

% save BLP and corrMaps
eval(['save ' date2 sess '_BLP BLP']);
eval(['save ' date2 sess '_corr Corr']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function franges = getFrequencyRanges()

franges = [ 1 4;...
            5 8;...
            9 14;...
           15 30;...
           35 80]';
           %300 499]';
% Delta [1-4] Theta [5-8] Alpha [9 14] Beta [15-30] Gamma [35 80] MUA [300 499]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%