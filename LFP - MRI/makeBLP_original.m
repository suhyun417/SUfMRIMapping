function BLP = makeBLP_original(LFP,fs)
% This function converts LFP input into band-limited power bands.
% USAGE:    BLP = makeBLP_original(LFP,Fs)
% INPUT:    LFP is the 2-dimensional LFP data, such that the first 
%           dimension represents time and the second dimension represents 
%           represents channel. fs is the (downsampled) sampling rate
%           resFile is the name of the resulting BLP structure
% OUTPUT:   BLP (structure with 3-dim matrices)
%           BLP.blp_dat w/ [channel x time x bandwidth]
%           BLP.fs_BLP
%           BLP.frange 


nyq = fs/2;
numchan = size(LFP,2);
frange = getFrequencyRanges;
nbw    = size(frange,2);

for chan = 1:numchan

    rdat = squeeze(LFP(:,chan));
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
    end

    BLP.blp_dat(chan,:,:)   = blp_dat; % channel x timepoint x bandwidth
    BLP.fs_BLP   = fs;
    BLP.frange   = frange;

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function franges = getFrequencyRanges()

franges = [ 2 6;...
            5 9;...
            8 14;...
           15 25;...
           25 40;...
           40 80]';
           %300 499]';
% Delta [1-4] Theta [5-8] Alpha [9 14] Beta [15-30] Gamma [35 80] MUA [300 499]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%