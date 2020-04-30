function dat =  build_blpDB(dat,band)
%
%  dat =  build_blpDB(dat,band)
%
% Inputs:
% Adds a .blp field to DAT. DAT is a Hanuman data structure that
% contains one LFP field. BAND is an N x 2 matrix. Each row of BAND
% is a 2-element vector specifying the upper and lower ends of a 
% frequency band in Hz. 
%
% Outputs:
% Dimensions of dat.blp is Ntrials x Ntimesamples x Nfreqencybands.
%
% Adds the following field to header:
% dat.h.blpfreq = band
%
% Example:
% band = [9 14; 15 30; 30 50; 50 100];
% dat = build_blpDB(dat,band);
% 
%
% 2014-jan-18
% Added capability for handling LFP inputs that contain NaNs. Currently
% assumes that all NaNs are consequtive and occur at the right-hand edge of
% LFP matrix.
%
% last modified 2014-jan-18
% dbtm

% get band-limited power    
Fs = 1000;
nq = Fs/2;
Nbands = size(band,1);
Ntrials = size(dat.c,1);
for b=1:Nbands
    [lb la] = cheby1(2,0.8,[band(b,:)]./nq);
    for t=1:Ntrials
        lfp = dat.lfp(t,:);
        mask = find(isnan(lfp));
        if ~isempty(mask)
            warning('LFP field contains NaNs. Assuming all NaNs at end...')
            if mask(end)~=size(dat.lfp,2)
                error('Ok, they weren''t all at the end.');
            end
            lfp(mask) = [];
            blpTrial = filtfilt(lb,la,lfp);
            blpTrial(mask) = NaN;
            blp(t,:,b) = blpTrial;
        else        
            blp(t,:,b)  = filtfilt(lb,la,dat.lfp(t,:));
        end
    end
end
blp = abs(blp);

dat.blp = blp;
dat.h.blpfreq = band;


% % smooth with adjustable Gaussian kernel
%sigma = 2*1000/band(2);
%mirror = fliplr(blp);
%padded = [mirror blp mirror];
%kwidth = -sigma*3:3*sigma;
%kernel = normpdf(kwidth,0,sigma);
%s = conv(padded,kernel);
%s(1:length(mirror)) = [];
%s(length(s)-length(mirror)+1:length(s)) = [];
%s(1:floor(length(kernel)/2)) = [];
%s(length(blp)+1:length(s)) = [];
%
%blp = s;
