function makeBLP
% This function converts LFP input into band-limited power bands.
% USAGE:    makeBLP
% INPUT:    LFP is the 3-dimensional LFP data, such that the first 
%           dimension represents time, the second dimension represents 
%           epoch (EPI volume), and the third dimension represents channel.
%           fs is the sampling rate; you should have your data already 
%           downsampled.
%           resFile is the name of the resulting BLP structure
% OUTPUT:   BLP (structure with 4-dim matrices)
%           BLP.blp_dat w/ [channel|time|epoch|bandwidth]
%           BLP.fs_BLP
%           BLP.frange

date = '11-08-08';  date2 = '110808';  sess = '_0002';  monkey = 'Varia';  fs = 250;

cd(['/einstein0/USRlab/projects/scholvinckm/data/' monkey '/inside scanner/' date '/Matlab']);
load([date2 sess '_LFP_epochs_new']);

LFP = permute(LFP,[1 3 2]);  % time x channel x epoch
nyq = fs/2;
numchan = size(LFP,2);
frange = getFrequencyRanges;
nbw    = size(frange,2);

for chan = 1:numchan

   rdat = squeeze(LFP(:,chan,:));
   filt_dat = zeros(size(rdat,1),size(rdat,2),nbw);
   blp_dat = zeros(size(rdat,1),size(rdat,2),nbw);
       
   for bw = 1:nbw
        fr                    = frange(:,bw);
        hpc                   = fr(1);                      % Hz
        lpc                   = fr(2);                      % Hz
        hWn                   = hpc/nyq;                    % between 0.0 and 1.0
        lWn                   = lpc/nyq;                    % between 0.0 and 1.0
        [lb,la]               = cheby1(2,0.8,[hWn lWn]);    % chebyshev bandpass
        fdat                  = filtfilt(lb,la,rdat);       % bandpass
        filt_dat(:,:,bw)      = fdat;
        blp_dat(:,:,bw)       = abs(fdat);                  % "power"
   end

%     % to diagnose problems
%     if 0
% 	  figure(3003);
% 	  clf
% 	  hold on
% 	  epoch = 1;
% 	  blptmp = squeeze(blp_dat(:,epoch,:));
% 	  filtered = squeeze(filt_dat(:,epoch,:));
% 	  size(repmat(mean(filtered,1),[size(filtered,1) 1]))
% 	  %plot(filtered-repmat(mean(filtered,1),[size(filtered,1) 1]))
% 	  plot(blptmp)
% 	  unfiltered = squeeze(rdat(:,epoch));	  
% 	  plot(unfiltered-mean(unfiltered),'k','LineWidth',2)
% 	  pause
%     end
	
    BLP.blp_dat(chan,:,:,:)   = blp_dat; % channel x timepoint x epoch x bandwidth
    BLP.fs_BLP   = fs;
    BLP.frange   = frange;
    
end

eval(['save ' date2 sess '_BLP_new BLP remove']);

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