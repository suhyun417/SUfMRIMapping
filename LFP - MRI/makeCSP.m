function makeCSP
% This function mean corrects LFP input, then computes the CSD from this 
% LFP and converts it into rectified band-limited power bands. Note: the 
% CSD calculation does not tolerate bad channels!
% USAGE:    makeCSP
% INPUT:    LFP is the 3-dimensional LFP data, such that the first 
%           dimension represents time, the second dimension represents 
%           epoch (EPI volume), and the third dimension represents channel.
%           fs is the sampling rate; you should have your data already 
%           downsampled.
%           resFile is the name of the resulting CSP structure
% OUTPUT:   CSP (structure with 4-dim matrices)
%           CSP.csp_dat w/ [channel|time|epoch|bandwidth]
%           CSP.fs_CSP
%           CSP.frange

date = '16-04-08';  date2 = '160408';  sess = '_0004';
monkey = 'Varia';
fs = 250;
el_pos = [0.1:0.2:3.2];      % electrode contact spacing
nyq = fs/2;
frange = getFrequencyRanges;
nbw    = size(frange,2);

cd(['/einstein0/USRlab/projects/scholvinckm/data/' monkey '/inside scanner/' date '/Matlab']);
load([date2 sess '_LFP_filter_epochs']);


% % mean correct LFP by subtracting overall mean of all epochs
% LFP = permute(LFP,[2 1 3]);  % epoch x time x channel
% max_uvolts = 150;
% for chan = 1:size(LFP,3)
%     avLFP = repmat(mean(mean(squeeze(LFP(:,:,chan)))),size(LFP,1),size(LFP,2)); 
%     corrLFP(:,:,chan) = squeeze(LFP(:,:,chan)) - avLFP;
%     if sum(abs(corrLFP(:,:,chan)')>max_uvolts), continue; end
% end
corrLFP = permute(LFP,[2 1 3]);   % epoch x time x channel

% COMPUTE CSD
[ntr,npt,nch] = size(corrLFP);
CSD = zeros(nch-2,npt,ntr);
tottr = size(corrLFP,1);

for trial = 1:tottr
    pot  = squeeze(corrLFP(trial,:,:))';
    cond     = 0.35;        
    cond_top = 0.35;        
    dt = 1/fs*1000;                 
    ddia = 0.5;            
    diam = ddia*1e-3;       
    N = length(el_pos);     
    d = mean(diff(el_pos));
    for i=1:N-2
        for j=1:N
            if (i == j-1)
                out(i,j) = -2/d^2;
            elseif (abs(i-j+1) == 1)
                out(i,j) = 1/d^2;
            else
                out(i,j) = 0;
            end;
        end;
    end;
    CSD(:,:,trial)  = -cond*out*pot; 
end


% BAND PASS FILTER AND RECTIFY CSD
CSD = permute(CSD,[2 1 3]);     % time x channel x epoch

for chan = 1:size(CSD,2)

   rdat = squeeze(CSD(:,chan,:));
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

    CSP.blp_dat(chan,:,:,:)   = blp_dat; % channel x timepoint x epoch x bandwidth
    CSP.fs_CSP   = fs;
    CSP.frange   = frange;
    
end

eval(['save ' date2 sess '_CSP_filter CSP CSD LFP remove']);


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