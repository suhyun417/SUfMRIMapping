function [sdf time] = getSdf(dat,win,varargin)
%
% [sdf time] = getSdf(dat,win)
% 
% Returns mean firing rate plot for trials in dat and corresponding time
% base.
% 
% sdf = getSdf(dat,win,sigma)
% 
% SIGMA Specifies smoothing kernel in ms (default = 100).
%
% Note that win and sigma should be in units consistent with dat.h. units.
% (either ms or sec).
%
% last modified 2013-apr-12
% dbtm


if ~isempty(varargin)
    %sp = varargin{1};
    sigma = varargin{1};
else
    %sp = 1;
    sigma = 100;
end

origUnits = dat.h.units;

if isequal(origUnits,'sec')
   %temporarily convert to ms
   dat = movieTimeScale(dat,'ms');
   win = win*1000;
   sigma = sigma*1000;
end

Ntrials = length(dat.s);
time = win(1):win(2);
spikes = zeros(Ntrials,length(time));

for t=1:Ntrials
    ts = round(dat.s{t});
    ts = ts(find((ts>win(1)) & ts<win(2)));
    ts = ts-win(1);
   spikes(t,ts) = 1000;
end
psth = mean(spikes,1);

% pad psth w/mirror-image of itself (to reduce edge-artifacts).
mirror = fliplr(psth);
padded = [mirror psth mirror];
% convolve psth with kernel
kwidth = -sigma*3:3*sigma;
kernel = openGaussian(kwidth,0,sigma);
s = conv(padded,kernel);
s(1:length(mirror)) = [];
s(length(s)-length(mirror)+1:length(s)) = [];
s(1:floor(length(kernel)/2)) = [];
s(length(psth)+1:length(s)) = [];
sdf = s;


if isequal(origUnits,'sec')
   % return to original condition.
   dat = movieTimeScale(dat,'sec');
   time = time/1000;
end













